#!.pixi/envs/default/bin/python

import duckdb
from pyfaidx import Fasta
from typing import Dict, List, Optional, Union
import pandas as pd
import logging
import os
from pathlib import Path
import threading
import time
from collections import OrderedDict
from Bio.SeqUtils import gc_fraction

from .fasta_utils import assemble_genome, get_genome_stats

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Maximum number of host FASTA files to keep open simultaneously.
# Prevents "too many open files" OS errors when using host-mapping mode.
MAX_HOST_FASTA_CACHE_SIZE = 100


def _fasta_key_function(header: str) -> str:
    """Extract the accession ID (first whitespace-delimited token) from a FASTA header."""
    parts = header.split()
    return parts[0] if parts else header


def _normalize_source_type(source_type: Optional[str]) -> str:
    """Normalize source type to canonical public/private labels."""
    normalized = str(source_type).strip().lower() if source_type is not None else ""
    return "private" if normalized == "private" else "public"


def _normalize_source_db(source_db: Optional[object]) -> Optional[str]:
    """Normalize source DB labels and discard empty/NaN-like values."""
    if source_db is None:
        return None
    normalized = str(source_db).strip()
    if not normalized or normalized.lower() == "nan":
        return None
    return normalized


def _should_rebuild_fai(fasta_path: Union[str, Path]) -> bool:
    """Return True when FASTA index is missing or older than the FASTA file."""
    fasta_path = Path(fasta_path)
    index_path = Path(str(fasta_path) + '.fai')
    if not index_path.exists():
        return True
    try:
        return index_path.stat().st_mtime < fasta_path.stat().st_mtime
    except OSError:
        return True


def _load_protein_fasta(path: str) -> "Fasta":
    """
    Load a protein FASTA file using full headers as keys.

    Protein FASTA headers have the phage accession as the first token
    (e.g. ">AE002163.1 CDS_1 hypothetical protein"), so multiple proteins
    from the same phage share the same first token.  Using only the first
    token as the key therefore causes pyfaidx to raise ``Duplicate key``.

    Using ``split_char='\\x00'`` (a character that never appears in FASTA
    headers) together with ``read_long_names=True`` makes pyfaidx read the
    *full* header line directly from the FASTA file and use it as the key,
    giving every protein sequence a unique, unambiguous identifier.
    """
    try:
        return Fasta(
            path,
            read_long_names=True,
            split_char='\x00',
        )
    except ValueError as e:
        if 'Duplicate key' in str(e):
            # The existing .fai was built with first-token keys; rebuild it.
            logging.warning(
                f"⚠️  Duplicate keys in protein FASTA index, rebuilding: {path}"
            )
            return Fasta(
                path,
                read_long_names=True,
                split_char='\x00',
                rebuild=True,
            )
        raise


def parse_where_clause(where_clause: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    """
    Parse a where_clause that may contain WHERE conditions, LIMIT, and/or OFFSET.
    
    Args:
        where_clause: The clause to parse (e.g., "LIMIT 100", "p.Length > 1000 LIMIT 50", 
                     "LIMIT 1000 OFFSET 5000", etc.)
    
    Returns:
        Tuple of (where_conditions, limit_offset_clause)
        - where_conditions: The WHERE conditions only (without LIMIT/OFFSET), or None
        - limit_offset_clause: The LIMIT/OFFSET clause only, or None
    
    Examples:
        >>> parse_where_clause("LIMIT 100")
        (None, "LIMIT 100")
        >>> parse_where_clause("p.Length > 1000 LIMIT 50")
        ("p.Length > 1000", "LIMIT 50")
        >>> parse_where_clause("LIMIT 1000 OFFSET 5000")
        (None, "LIMIT 1000 OFFSET 5000")
        >>> parse_where_clause("p.GC > 0.5")
        ("p.GC > 0.5", None)
    """
    if not where_clause:
        return None, None
    
    # Normalize whitespace
    clause = ' '.join(where_clause.split())
    
    # Case-insensitive search for LIMIT keyword
    clause_upper = clause.upper()
    limit_pos = clause_upper.find(' LIMIT ')
    
    # Handle case where LIMIT is at the start
    if clause_upper.startswith('LIMIT '):
        limit_pos = 0
    
    if limit_pos == -1:
        # No LIMIT clause found
        return clause.strip(), None
    elif limit_pos == 0:
        # LIMIT is at the start, no WHERE conditions
        return None, clause.strip()
    else:
        # Split at LIMIT position
        where_part = clause[:limit_pos].strip()
        limit_part = clause[limit_pos:].strip()
        
        return where_part if where_part else None, limit_part if limit_part else None


class SequenceRetriever:
    """
    Retrieve sequences from indexed FASTA files based on DuckDB queries
    
    Features:
    - Lazy loading with background thread support
    - Query-based sequence retrieval from DuckDB
    - Direct ID-based retrieval
    - Batch processing support
    - Memory-efficient streaming
    """
    
    def __init__(self, db_path: str, phage_fasta_path: str, protein_fasta_path: str, 
                 host_fasta_path: Optional[str] = None, host_mapping_path: Optional[str] = None,
                 private_phage_mapping_path: Optional[str] = None,
                 preload: bool = True):
        """
        Initialize SequenceRetriever with lazy FASTA loading
        
        Args:
            db_path: Path to DuckDB database
            phage_fasta_path: Path to indexed phage FASTA file (public sequences only)
            protein_fasta_path: Path to indexed protein FASTA file
            host_fasta_path: Path to indexed host FASTA file (DEPRECATED - use host_mapping_path)
            host_mapping_path: Path to JSON mapping file for individual host FASTA files
            private_phage_mapping_path: Path to JSON mapping ``source_db → phage.fasta`` for
                private datasets.  When provided, SequenceRetriever routes private-phage
                sequence lookups to the matching per-source FASTA instead of the public
                ``phage_fasta_path``.  Pass the value of ``config["private_phage_mapping"]``
                from the Snakemake config.
            preload: If True, load FASTA files in background thread (default: True)
        """
        # Validate paths
        if not Path(db_path).exists():
            raise FileNotFoundError(f"Database not found: {db_path}")
        if not Path(phage_fasta_path).exists():
            raise FileNotFoundError(f"Phage FASTA not found: {phage_fasta_path}")
        if not Path(protein_fasta_path).exists():
            raise FileNotFoundError(f"Protein FASTA not found: {protein_fasta_path}")
        
        # Check for index files
        phage_index = Path(str(phage_fasta_path) + '.fai')
        protein_index = Path(str(protein_fasta_path) + '.fai')
        
        logging.info(f"📂 Checking FASTA index files:")
        logging.info(f"   Phage index: {phage_index.exists()} ({phage_index.stat().st_size / 1024:.1f} KB)")
        logging.info(f"   Protein index: {protein_index.exists()} ({protein_index.stat().st_size / 1024:.1f} KB)")
        
        if not phage_index.exists():
            raise FileNotFoundError(f"Phage FASTA index not found: {phage_index}")
        if not protein_index.exists():
            raise FileNotFoundError(f"Protein FASTA index not found: {protein_index}")
        
        # Initialize private phage mapping (source_db → per-source phage.fasta path).
        # Private phages are NOT in all_phages.fasta; they live in per-source FASTA
        # files under /data/intermediate/fasta/private/phages/<source_db>/phage.fasta.
        import json as _json
        self._private_phage_mapping: Optional[Dict[str, str]] = None
        self._private_phage_source_lookup: Dict[str, str] = {}
        self._private_phage_fasta_cache: OrderedDict = OrderedDict()
        self._private_phage_lock = threading.Lock()

        if private_phage_mapping_path:
            _ppm = Path(private_phage_mapping_path)
            if _ppm.exists():
                with _ppm.open("r") as _f:
                    self._private_phage_mapping = _json.load(_f)
                for _source_key in self._private_phage_mapping.keys():
                    _normalized_source = _normalize_source_db(_source_key)
                    if _normalized_source:
                        _norm_key = _normalized_source.lower()
                        _existing = self._private_phage_source_lookup.get(_norm_key)
                        if _existing is not None and _existing != _source_key:
                            logging.warning(
                                "⚠️ Duplicate normalized private source key '%s' for '%s' and '%s'; using '%s'",
                                _norm_key,
                                _existing,
                                _source_key,
                                _existing,
                            )
                        else:
                            self._private_phage_source_lookup[_norm_key] = _source_key
                logging.info(
                    f"📂 Loaded private phage mapping for "
                    f"{len(self._private_phage_mapping)} sources: {list(self._private_phage_mapping.keys())}"
                )
            else:
                logging.debug(f"Private phage mapping not found: {private_phage_mapping_path}")

        # Initialize host data handling
        self._host_fasta_path = host_fasta_path  # Legacy single-file mode
        self._host_mapping_path = host_mapping_path  # New mapping mode
        self._host_mapping = None  # Mapping from Host_ID to file path
        self._host_fasta_cache = OrderedDict()  # LRU cache of loaded Fasta objects per host
        self._host_fasta = None  # For legacy single-file mode
        self._host_lock = threading.Lock()
        self._host_count = None
        self._has_host_data = False
        self._use_host_mapping = False
        
        # Check if using new mapping-based approach
        if host_mapping_path:
            if Path(host_mapping_path).exists():
                logging.info(f"📂 Using host mapping file: {host_mapping_path}")
                self._has_host_data = True
                self._use_host_mapping = True
                # Load mapping file
                import json
                with open(host_mapping_path, 'r') as f:
                    self._host_mapping = json.load(f)
                self._host_count = len(self._host_mapping)
                logging.info(f"   Loaded mapping for {self._host_count} hosts")
            else:
                logging.warning(f"⚠️  Host mapping file not found: {host_mapping_path}")
        # Fallback to legacy single-file mode
        elif host_fasta_path:
            if Path(host_fasta_path).exists():
                host_index = Path(str(host_fasta_path) + '.fai')
                if host_index.exists():
                    logging.info(f"   Host index: {host_index.exists()} ({host_index.stat().st_size / 1024:.1f} KB)")
                    self._has_host_data = True
                else:
                    logging.warning(f"⚠️  Host FASTA index not found: {host_index}")
            else:
                logging.warning(f"⚠️  Host FASTA not found: {host_fasta_path}")
        
        # Initialize database connection (fast)
        logging.info(f"📂 Connecting to database: {db_path}")
        self.conn = duckdb.connect(db_path, read_only=True)
        
        # Store paths for lazy loading
        self._phage_fasta_path = phage_fasta_path
        self._protein_fasta_path = protein_fasta_path
        
        # Lazy-loaded FASTA objects
        self._phage_fasta = None
        self._protein_fasta = None
        
        # Thread synchronization
        self._phage_lock = threading.Lock()
        self._protein_lock = threading.Lock()
        self._loading_complete = threading.Event()
        
        # Stats
        self._phage_count = None
        self._protein_count = None
        
        if preload:
            # Start background loading
            logging.info("🔄 Starting background FASTA loading...")
            self._load_thread = threading.Thread(target=self._preload_fasta, daemon=True)
            self._load_thread.start()
            logging.info("✅ Initialization complete (FASTA loading in background)")
        else:
            logging.info("✅ Initialization complete (FASTA files will load on first use)")
    
    def _preload_fasta(self):
        """Background task to load FASTA files"""
        start_total = time.time()
        
        try:
            # Load phage FASTA
            logging.info(f"🔄 [Background] Loading phage FASTA: {self._phage_fasta_path}")
            start = time.time()
            
            with self._phage_lock:
                self._phage_fasta = Fasta(
                    self._phage_fasta_path,
                    rebuild=False,
                    key_function=_fasta_key_function
                )
                self._phage_count = len(self._phage_fasta.keys())
            
            elapsed = time.time() - start
            logging.info(f"   ✅ Phage FASTA loaded in {elapsed:.2f}s ({self._phage_count:,} sequences)")
            
            # Load protein FASTA
            logging.info(f"🔄 [Background] Loading protein FASTA: {self._protein_fasta_path}")
            start = time.time()
            
            with self._protein_lock:
                self._protein_fasta = _load_protein_fasta(self._protein_fasta_path)
                self._protein_count = len(self._protein_fasta.keys())
            
            elapsed = time.time() - start
            logging.info(f"   ✅ Protein FASTA loaded in {elapsed:.2f}s ({self._protein_count:,} sequences)")
            
            # Load host FASTA if available (only for legacy single-file mode)
            if self._has_host_data and self._host_fasta_path and not self._use_host_mapping:
                logging.info(f"🔄 [Background] Loading host FASTA: {self._host_fasta_path}")
                start = time.time()
                
                with self._host_lock:
                    self._host_fasta = Fasta(
                        self._host_fasta_path,
                        rebuild=False,
                        key_function=_fasta_key_function
                    )
                    self._host_count = len(self._host_fasta.keys())
                
                elapsed = time.time() - start
                logging.info(f"   ✅ Host FASTA loaded in {elapsed:.2f}s ({self._host_count:,} sequences)")
            elif self._use_host_mapping:
                logging.info(f"   ℹ️  Using on-demand loading for {self._host_count:,} individual host files")
            
            # Mark as complete
            self._loading_complete.set()
            
            total_elapsed = time.time() - start_total
            logging.info(f"🎉 All FASTA files loaded in {total_elapsed:.2f}s")
            
            # Sample keys
            sample_phage = list(self._phage_fasta.keys())[:3]
            sample_protein = list(self._protein_fasta.keys())[:3]
            
            logging.info(f"🔍 Sample phage keys:")
            for key in sample_phage:
                logging.info(f"   - '{key[:80]}...'")
            
            logging.info(f"🔍 Sample protein keys:")
            for key in sample_protein:
                logging.info(f"   - '{key[:80]}...'")
            
            if self._host_fasta:
                sample_host = list(self._host_fasta.keys())[:3]
                logging.info(f"🔍 Sample host keys:")
                for key in sample_host:
                    logging.info(f"   - '{key[:80]}...'")
                
        except Exception as e:
            logging.error(f"❌ Error loading FASTA files: {e}")
            import traceback
            traceback.print_exc()
    
    @property
    def phage_fasta(self):
        """Get phage FASTA, loading if necessary"""
        if self._phage_fasta is None:
            with self._phage_lock:
                if self._phage_fasta is None:  # Double-check locking
                    logging.info(f"📂 Loading phage FASTA on-demand: {self._phage_fasta_path}")
                    start = time.time()
                    self._phage_fasta = Fasta(
                        self._phage_fasta_path,
                        rebuild=False,
                        key_function=_fasta_key_function
                    )
                    elapsed = time.time() - start
                    logging.info(f"   ✅ Loaded in {elapsed:.2f}s")
        return self._phage_fasta
    
    @property
    def protein_fasta(self):
        """Get protein FASTA, loading if necessary"""
        if self._protein_fasta is None:
            with self._protein_lock:
                if self._protein_fasta is None:  # Double-check locking
                    logging.info(f"📂 Loading protein FASTA on-demand: {self._protein_fasta_path}")
                    start = time.time()
                    self._protein_fasta = _load_protein_fasta(self._protein_fasta_path)
                    elapsed = time.time() - start
                    logging.info(f"   ✅ Loaded in {elapsed:.2f}s")
        return self._protein_fasta
    
    @property
    def host_fasta(self):
        """
        Get host FASTA, loading if necessary
        
        DEPRECATED: This property is maintained for backward compatibility
        but is not available when using host_mapping_path. Use 
        get_host_sequence() method instead for individual host access.
        """
        if not self._has_host_data:
            raise ValueError("Host FASTA not configured - pass host_fasta_path or host_mapping_path to __init__")
        
        # If using mapping mode, this operation is not supported
        if self._use_host_mapping:
            raise ValueError(
                "Direct access to host_fasta is not available when using host_mapping_path. "
                "The host genomes are stored as individual files for efficiency. "
                "Use get_host_sequence(host_id) method instead to load individual host files on-demand."
            )
        
        # Legacy single-file mode
        if self._host_fasta is None:
            with self._host_lock:
                if self._host_fasta is None:  # Double-check locking
                    logging.info(f"📂 Loading host FASTA on-demand: {self._host_fasta_path}")
                    start = time.time()
                    self._host_fasta = Fasta(
                        self._host_fasta_path,
                        rebuild=False,
                        key_function=_fasta_key_function
                    )
                    elapsed = time.time() - start
                    logging.info(f"   ✅ Loaded in {elapsed:.2f}s")
        return self._host_fasta
    
    def _get_host_fasta_for_id(self, host_id: str) -> Fasta:
        """
        Get Fasta object for a specific host ID (used in mapping mode)
        
        Args:
            host_id: Host identifier
            
        Returns:
            Fasta object for the host genome
            
        Raises:
            KeyError: If host_id not found in mapping
            FileNotFoundError: If host file doesn't exist
        """
        if not self._use_host_mapping:
            # Should not be called in legacy mode
            raise RuntimeError("This method is only for host mapping mode")
        
        # Get file path from mapping (check before acquiring lock)
        if host_id not in self._host_mapping:
            raise KeyError(f"Host ID '{host_id}' not found in mapping")
        
        fasta_path = self._host_mapping[host_id]
        
        # Load the fasta file with LRU cache management
        with self._host_lock:
            if host_id in self._host_fasta_cache:
                # Move to end (most recently used)
                self._host_fasta_cache.move_to_end(host_id)
                return self._host_fasta_cache[host_id]

            fasta_path = self._resolve_host_fasta_path(host_id, fasta_path)
            
            # Evict oldest entry if cache is full
            if len(self._host_fasta_cache) >= MAX_HOST_FASTA_CACHE_SIZE:
                oldest_id, oldest_fasta = self._host_fasta_cache.popitem(last=False)
                if hasattr(oldest_fasta, 'close'):
                    try:
                        oldest_fasta.close()
                    except Exception as e:
                        logging.debug(f"Error closing evicted host FASTA for {oldest_id}: {e}")
            
            rebuild = _should_rebuild_fai(fasta_path)
            if rebuild:
                logging.info(f"Creating index for {fasta_path}")
            
            logging.debug(f"Loading host FASTA for {host_id}: {fasta_path}")
            fasta_obj = Fasta(
                fasta_path,
                rebuild=rebuild,
                key_function=_fasta_key_function
            )
            self._host_fasta_cache[host_id] = fasta_obj
        
        return self._host_fasta_cache[host_id]

    def _resolve_host_fasta_path(self, host_id: str, mapped_path: str) -> str:
        """
        Resolve host FASTA path from mapping, with fallback search for stale paths.

        If the mapped path no longer exists (e.g. the private-data mount point changed),
        this searches ``<root>/*/hosts/<filename>`` under known private-data roots.
        """
        path = Path(mapped_path)
        if path.exists():
            return str(path)

        # Handle relative mapping entries against CWD.
        if not path.is_absolute():
            cwd_candidate = Path.cwd() / path
            if cwd_candidate.exists():
                return str(cwd_candidate)

        # Build fallback roots where private source folders may be available.
        candidate_roots = []
        env_private_root = os.getenv("PBI_PRIVATE_DATA_DIR")
        if env_private_root:
            candidate_roots.append(Path(env_private_root))
        candidate_roots.extend([
            Path("/private-data"),
            Path(__file__).resolve().parents[2] / "private_data",
            Path.cwd() / "private_data",
        ])

        expected_source = None
        if len(path.parts) >= 3 and path.parent.name == "hosts":
            expected_source = path.parent.parent.name or None

        seen = set()
        for root in candidate_roots:
            root_str = str(root.expanduser().resolve(strict=False))
            if root_str in seen:
                continue
            seen.add(root_str)
            if not root.exists():
                continue
            if root.is_file():
                continue

            matches = []

            # Prefer the source folder from the stale path to avoid accidental picks
            # from hidden cache folders (e.g. .pbi) when both contain same filename.
            if expected_source:
                preferred = root / expected_source / "hosts" / path.name
                if preferred.exists():
                    matches.append(preferred)

            if not matches:
                # Then scan non-hidden source folders.
                for source_dir in sorted(
                    (p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".")),
                    key=lambda p: p.name,
                ):
                    # Deterministic ordering keeps resolution stable across runs.
                    candidate = source_dir / "hosts" / path.name
                    if candidate.exists():
                        matches.append(candidate)

            if not matches:
                # Last-resort compatibility path: include hidden directories.
                matches = sorted(root.glob(f"*/hosts/{path.name}"))

            if matches:
                resolved = matches[0]
                logging.warning(
                    "⚠️ Resolved missing host mapping path for %s: %s -> %s",
                    host_id,
                    mapped_path,
                    resolved,
                )
                self._host_mapping[host_id] = str(resolved)
                return str(resolved)

        return mapped_path

    def _resolve_private_phage_fasta_path(self, source_db: str, mapped_path: str) -> str:
        """
        Resolve private phage FASTA path from mapping, with fallback for stale mounts.

        If the mapped path no longer exists (e.g. ``/private-data`` moved to another
        mount path), this searches common private-data roots for a matching
        ``<source_db>/phage.fasta`` and updates the in-memory mapping.

        Returns:
            Resolved existing FASTA path when found; otherwise returns ``mapped_path``
            unchanged.
        """
        path = Path(mapped_path)
        if path.exists():
            return str(path)

        # Handle relative mapping entries against CWD.
        if not path.is_absolute():
            cwd_candidate = Path.cwd() / path
            if cwd_candidate.exists():
                return str(cwd_candidate)

        candidate_roots = []
        env_private_root = os.getenv("PBI_PRIVATE_DATA_DIR")
        if env_private_root:
            candidate_roots.append(Path(env_private_root))
        candidate_roots.extend(
            [
                Path("/private-data"),
                Path(__file__).resolve().parents[2] / "private_data",
                Path.cwd() / "private_data",
            ]
        )

        seen = set()
        for root in candidate_roots:
            root_str = str(root.expanduser().absolute())
            if root_str in seen:
                continue
            seen.add(root_str)
            if not root.exists() or not root.is_dir():
                continue

            matches = []

            # Preserve path suffix when remapping a stale /private-data mount.
            mapped_norm = path.as_posix()
            private_data_prefix = "/private-data/"
            if private_data_prefix in mapped_norm:
                remapped_suffix = mapped_norm.split(private_data_prefix, 1)[1]
                remapped = root / Path(remapped_suffix)
                if remapped.exists():
                    matches.append(remapped)

            # Common layout for source-mounted private datasets.
            preferred = root / source_db / "phage.fasta"
            if preferred.exists():
                matches.append(preferred)

            # Common layout for copied/indexed per-source private phage FASTAs.
            intermediate = root / "private_phage_genomes_intermediate" / source_db / "phage.fasta"
            if intermediate.exists():
                matches.append(intermediate)

            if not matches:
                # Last resort: shallow lookup for */<source_db>/phage.fasta.
                # Keep sorted order so fallback resolution is deterministic.
                for candidate in sorted(root.glob(f"*/{source_db}/phage.fasta")):
                    try:
                        relative_parts = candidate.relative_to(root).parts
                    except ValueError:
                        continue
                    if any(part.startswith(".") for part in relative_parts):
                        continue
                    matches.append(candidate)
                    break

            if matches:
                resolved = matches[0]
                logging.warning(
                    "⚠️ Resolved missing private phage mapping path for source %s: %s -> %s",
                    source_db,
                    mapped_path,
                    resolved,
                )
                self._private_phage_mapping[source_db] = str(resolved)
                return str(resolved)

        return mapped_path

    def _resolve_private_source_db(self, source_db: Optional[object]) -> Optional[str]:
        """Resolve a DB-provided source label to a canonical private mapping key."""
        normalized_source = _normalize_source_db(source_db)
        if normalized_source is None or not self._private_phage_mapping:
            return None
        if normalized_source in self._private_phage_mapping:
            return normalized_source
        return self._private_phage_source_lookup.get(normalized_source.lower())

    def _get_private_phage_fasta(self, source_db: str) -> "Fasta":
        """
        Load (and LRU-cache) the phage FASTA for a private source.

        Private phage FASTAs live in per-source directories under
        ``private_phage_genomes_intermediate/<source_db>/phage.fasta``.
        The mapping from ``source_db`` to the concrete FASTA path was built by
        ``prepare_private_sequence_artifacts`` and is loaded at construction time.

        Args:
            source_db: Source database identifier (e.g. ``"test_private"``).

        Returns:
            ``pyfaidx.Fasta`` object keyed by first-token of the FASTA header
            (same convention as the public phage FASTA).

        Raises:
            KeyError: If *source_db* is not in the private phage mapping.
            FileNotFoundError: If the FASTA file no longer exists on disk.
        """
        resolved_source_db = self._resolve_private_source_db(source_db)
        if (
            self._private_phage_mapping is None
            or resolved_source_db is None
            or resolved_source_db not in self._private_phage_mapping
        ):
            raise KeyError(f"Private phage source '{source_db}' not found in private_phage_mapping")

        with self._private_phage_lock:
            if resolved_source_db in self._private_phage_fasta_cache:
                self._private_phage_fasta_cache.move_to_end(resolved_source_db)
                return self._private_phage_fasta_cache[resolved_source_db]

            # Evict oldest entry when the cache is full (same limit as host cache)
            if len(self._private_phage_fasta_cache) >= MAX_HOST_FASTA_CACHE_SIZE:
                oldest_id, oldest_fasta = self._private_phage_fasta_cache.popitem(last=False)
                if hasattr(oldest_fasta, "close"):
                    try:
                        oldest_fasta.close()
                    except Exception:
                        pass

            fasta_path = self._private_phage_mapping[resolved_source_db]
            fasta_path = self._resolve_private_phage_fasta_path(resolved_source_db, fasta_path)
            rebuild = _should_rebuild_fai(fasta_path)
            if rebuild:
                logging.info("Creating index for private phage FASTA: %s", fasta_path)

            fasta_obj = Fasta(fasta_path, rebuild=rebuild, key_function=_fasta_key_function)
            self._private_phage_fasta_cache[resolved_source_db] = fasta_obj

        return self._private_phage_fasta_cache[resolved_source_db]

    def _get_phage_sequence(
        self, phage_id: str, source_db: Optional[str] = None, source_type: Optional[str] = None
    ) -> Optional[str]:
        """
        Retrieve a single phage sequence, routing private phages to their source FASTA.

        If *source_type* is ``"private"`` and a private phage mapping is configured,
        the sequence is fetched from the per-source ``phage.fasta``.  Otherwise the
        public ``all_phages.fasta`` is used.

        Args:
            phage_id: Phage identifier (first token of the FASTA header).
            source_db: Source database name.  Required when routing private phages.
            source_type: ``"private"`` or ``"public"`` (or ``None`` for unknown).

        Returns:
            Sequence string, or ``None`` if not found.
        """
        normalized_source_type = _normalize_source_type(source_type)
        resolved_source_db = self._resolve_private_source_db(source_db)
        # Mapping membership is treated as authoritative for private routing even if
        # source_type is missing/mislabelled in upstream metadata.
        is_private = normalized_source_type == "private" or resolved_source_db is not None

        if is_private and self._private_phage_mapping and resolved_source_db:
            try:
                fasta_obj = self._get_private_phage_fasta(resolved_source_db)
                return str(fasta_obj[phage_id][:].seq)
            except (KeyError, Exception) as exc:
                logging.debug(
                    "Private phage %s not found in source %s: %s",
                    phage_id,
                    resolved_source_db,
                    exc,
                )

        try:
            return str(self.phage_fasta[phage_id][:].seq)
        except KeyError:
            # If private mapping is present but source_db was missing/misaligned, scan all
            # private sources as a last-resort fallback.
            if self._private_phage_mapping:
                for sdb in self._private_phage_mapping:
                    try:
                        fasta_obj = self._get_private_phage_fasta(sdb)
                        if phage_id in fasta_obj:
                            logging.debug(
                                "Found private phage %s in fallback source %s", phage_id, sdb
                            )
                            return str(fasta_obj[phage_id][:].seq)
                    except Exception:
                        pass
            return None

    def get_host_sequence(self, host_id: str, contig_mode: str = "first") -> str:
        """
        Get sequence for a specific host ID.

        Works with both legacy single-file mode and new mapping mode.
        In mapping mode, loads individual host files on-demand.

        Args:
            host_id: Host identifier.
            contig_mode: How to handle multi-contig FASTA files.

                * ``"first"`` *(default)* – return only the first (longest)
                  contig.  Preserves the original single-string behaviour.
                * ``"concat"`` – concatenate all contigs into one string
                  (contigs sorted by length desc, then header asc).  Useful
                  when the host genome is fragmented across scaffolds.

                .. note::
                    For the legacy single-merged-file mode the host entry
                    is always a single record, so both modes are equivalent.

        Returns:
            Host sequence as a single string.

        Raises:
            KeyError: If *host_id* is not found.
            ValueError: If *contig_mode* is not a recognised value.
        """
        if not self._has_host_data:
            raise ValueError("Host FASTA not configured")

        if contig_mode not in ("first", "concat"):
            raise ValueError(
                f"contig_mode must be 'first' or 'concat', got '{contig_mode}'."
            )

        if self._use_host_mapping:
            # New mapping mode - load individual file
            fasta_obj = self._get_host_fasta_for_id(host_id)
            if not fasta_obj.keys():
                raise KeyError(f"No sequences found in host file for {host_id}")
            return assemble_genome(fasta_obj, mode=contig_mode)
        else:
            # Legacy mode - use single merged file
            return str(self.host_fasta[host_id][:].seq)
    
    def get_host_genome(
        self,
        host_id: str,
        mode: str = "concat",
        gap: int = 0,
        order: str = "length_desc",
    ) -> Union[str, List[str], Dict[str, str]]:
        """Retrieve the full genome for a host, handling multi-contig FASTA files.

        This is the preferred method when a host genome is fragmented across
        multiple scaffolds or chromosomes.

        In mapping mode (*host_mapping_path* was provided at construction time),
        each host has its own FASTA file that may contain many contigs.  This
        method assembles them according to *mode* and *order*.

        In legacy single-file mode, the host entry is always a single record,
        so *mode* has no practical effect (all modes return equivalent results).

        Args:
            host_id: Host identifier (must exist in the host mapping).
            mode: Assembly mode – see :func:`~pbi.fasta_utils.assemble_genome`
                for full documentation.

                * ``"concat"`` *(default)* – all contigs joined into one string.
                * ``"first"`` – only the first contig (longest by default).
                * ``"list"`` – list of per-contig strings.
                * ``"dict"`` – ``{header: sequence}`` mapping.

            gap: Number of ``N`` characters to insert between contigs when
                *mode* is ``"concat"``.  Default is ``0``.
            order: Contig ordering.  ``"length_desc"`` *(default)* or
                ``"file"``.  See :func:`~pbi.fasta_utils.assemble_genome`.

        Returns:
            * ``str`` when *mode* is ``"concat"`` or ``"first"``.
            * ``list[str]`` when *mode* is ``"list"``.
            * ``dict[str, str]`` when *mode* is ``"dict"``.

        Raises:
            ValueError: If host data is not configured or arguments are invalid.
            KeyError: If *host_id* is not found.

        Example:
            >>> # Concatenated genome (default):
            >>> seq = retriever.get_host_genome("GCF_000005845")

            >>> # With 100-N gap between scaffolds:
            >>> seq = retriever.get_host_genome("GCF_000005845", gap=100)

            >>> # List of individual contig sequences:
            >>> contigs = retriever.get_host_genome("GCF_000005845", mode="list")

            >>> # Statistics only (no assembly):
            >>> stats = retriever.get_host_genome_stats("GCF_000005845")
        """
        if not self._has_host_data:
            raise ValueError("Host FASTA not configured")

        if self._use_host_mapping:
            fasta_obj = self._get_host_fasta_for_id(host_id)
            if not fasta_obj.keys():
                raise KeyError(f"No sequences found in host file for {host_id}")
            return assemble_genome(fasta_obj, mode=mode, gap=gap, order=order)
        else:
            # Legacy mode – always a single record; all modes are equivalent.
            seq = str(self.host_fasta[host_id][:].seq)
            if mode == "list":
                return [seq]
            if mode == "dict":
                return {host_id: seq}
            return seq

    def get_host_genome_stats(
        self,
        host_id: str,
        order: str = "length_desc",
    ) -> Dict[str, object]:
        """Return contig statistics for a host genome FASTA.

        Useful for inspecting assembly fragmentation without loading the full
        sequence data.

        Args:
            host_id: Host identifier.
            order: Contig ordering for the returned *lengths* list.
                ``"length_desc"`` *(default)* or ``"file"``.

        Returns:
            dict with keys:

            * ``"contig_count"`` (:class:`int`) – number of records.
            * ``"lengths"`` (:class:`list[int]`) – per-contig lengths.
            * ``"total_length"`` (:class:`int`) – sum of all lengths.

        Raises:
            ValueError: If host data is not configured.
            KeyError: If *host_id* is not found.

        Example:
            >>> stats = retriever.get_host_genome_stats("GCF_000005845")
            >>> print(stats["contig_count"])   # e.g. 7
            >>> print(stats["total_length"])   # e.g. 5400000
        """
        if not self._has_host_data:
            raise ValueError("Host FASTA not configured")

        if self._use_host_mapping:
            fasta_obj = self._get_host_fasta_for_id(host_id)
            return get_genome_stats(fasta_obj, order=order)
        else:
            seq = str(self.host_fasta[host_id][:].seq)
            length = len(seq)
            return {"contig_count": 1, "lengths": [length], "total_length": length}

    def get_phage_genome(
        self,
        phage_id: str,
        mode: str = "concat",
        gap: int = 0,
        order: str = "length_desc",
    ) -> Union[str, List[str], Dict[str, str]]:
        """Retrieve the full genome for a phage, handling multi-contig cases.

        Most phage genomes are single-contig; this method handles the rare
        cases where a phage FASTA entry is split across multiple records.

        When there is exactly one record for *phage_id*, all modes return
        equivalent results (a single sequence string / single-element
        list or dict).

        .. note::
            This method operates on the **per-phage record** looked up by
            *phage_id* key.  It does *not* attempt to group records that share
            a common prefix.  For the typical single-record phage case,
            ``get_phage_genome(id)`` is equivalent to ``get_phage_sequence(id)``.

        Args:
            phage_id: Phage identifier (must exist as a key in the phage FASTA).
            mode: Assembly mode.  See :func:`~pbi.fasta_utils.assemble_genome`.
                Default is ``"concat"``.
            gap: Gap N-characters between contigs for ``mode="concat"``.
                Default is ``0``.
            order: Contig ordering.  ``"length_desc"`` *(default)* or
                ``"file"``.

        Returns:
            * ``str`` when *mode* is ``"concat"`` or ``"first"``.
            * ``list[str]`` when *mode* is ``"list"``.
            * ``dict[str, str]`` when *mode* is ``"dict"``.

        Raises:
            KeyError: If *phage_id* is not found in the phage FASTA.

        Example:
            >>> seq = retriever.get_phage_genome("NC_000866")
            >>> # For a typical single-contig phage this is the same as:
            >>> seq = str(retriever.phage_fasta["NC_000866"][:].seq)
        """
        # Query the DB for source info so private phages are routed to their FASTA.
        source_db: Optional[str] = None
        source_type_val: Optional[str] = None
        if self._private_phage_mapping:
            try:
                row = self.conn.execute(
                    "SELECT Source_DB, source_type FROM fact_phages WHERE Phage_ID = ? LIMIT 1",
                    [phage_id],
                ).fetchone()
                if row:
                    source_db, source_type_val = row
            except Exception:
                pass

        seq = self._get_phage_sequence(phage_id, source_db=source_db, source_type=source_type_val)
        if seq is None:
            raise KeyError(f"Phage ID '{phage_id}' not found in phage FASTA.")

        # Phage FASTA is keyed by accession – one key → one record.
        # Wrap into a minimal single-entry container to honour the mode API.
        if mode == "list":
            return [seq]
        if mode == "dict":
            return {phage_id: seq}
        # "first" and "concat" both reduce to the single sequence.
        return seq

    def wait_until_ready(self, timeout: Optional[float] = None):
        """
        Wait for background loading to complete
        
        Args:
            timeout: Maximum seconds to wait (None = wait forever)
        
        Returns:
            bool: True if loading completed, False if timeout
        """
        if not hasattr(self, '_load_thread'):
            return True  # No background loading was started
        
        logging.info("⏳ Waiting for FASTA loading to complete...")
        result = self._loading_complete.wait(timeout=timeout)
        
        if result:
            logging.info("✅ FASTA loading complete")
        else:
            logging.warning(f"⚠️ Timeout after {timeout}s - FASTA may still be loading")
        
        return result
    
    def is_ready(self) -> bool:
        """Check if FASTA files are loaded"""
        return self._loading_complete.is_set()
    
    def get_phage_sequences(self, query: str, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get phage sequences based on SQL query
        
        Note: Will wait for phage FASTA to load if background loading is still in progress
        """
        # Ensure phage FASTA is loaded (property handles this)
        _ = self.phage_fasta  # Trigger loading if needed
        
        logging.info(f"🔍 Executing query: {query[:100]}...")
        
        if limit:
            query = f"{query} LIMIT {limit}"
        
        result = self.conn.execute(query).fetchdf()
        
        if 'Phage_ID' not in result.columns:
            raise ValueError("Query must return 'Phage_ID' column")
        
        phage_ids = result['Phage_ID'].tolist()
        logging.info(f"📊 Retrieved {len(phage_ids):,} Phage IDs from query")
        
        return self._fetch_phage_sequences(phage_ids)
    
    def get_protein_sequences(self, query: str, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get protein sequences based on SQL query
        
        Note: Will wait for protein FASTA to load if background loading is still in progress
        """
        # Ensure protein FASTA is loaded (property handles this)
        _ = self.protein_fasta  # Trigger loading if needed
        
        logging.info(f"🔍 Executing query: {query[:100]}...")
        
        if limit:
            query = f"{query} LIMIT {limit}"
        
        result = self.conn.execute(query).fetchdf()
        
        if 'Protein_ID' not in result.columns:
            raise ValueError("Query must return 'Protein_ID' column")
        
        protein_ids = result['Protein_ID'].tolist()
        logging.info(f"📊 Retrieved {len(protein_ids):,} Protein IDs from query")
        
        return self._fetch_protein_sequences(protein_ids)
    
    def get_stats(self) -> Dict:
        """
        Get database and FASTA statistics
        
        Note: Will wait for FASTA loading to complete
        """
        # Wait for background loading if still running
        if hasattr(self, '_load_thread') and not self.is_ready():
            self.wait_until_ready(timeout=300)  # 5 minute timeout
        
        # Ensure FASTA files are loaded
        _ = self.phage_fasta
        _ = self.protein_fasta
        
        stats = {
            'database': {
                'phages': self.conn.execute("SELECT COUNT(*) FROM fact_phages").fetchone()[0],
                'proteins': self.conn.execute("SELECT COUNT(*) FROM dim_proteins").fetchone()[0],
            },
            'fasta': {
                'phages': self._phage_count if self._phage_count else len(self.phage_fasta.keys()),
                'proteins': self._protein_count if self._protein_count else len(self.protein_fasta.keys()),
            }
        }
        
        # Add host stats if available
        if self._has_host_data:
            try:
                stats['database']['hosts'] = self.conn.execute("SELECT COUNT(*) FROM dim_hosts").fetchone()[0]
                stats['database']['phage_host_associations'] = self.conn.execute(
                    "SELECT COUNT(*) FROM phage_host_associations"
                ).fetchone()[0]
                
                # For mapping mode, use the count from the mapping
                if self._use_host_mapping:
                    stats['fasta']['hosts'] = self._host_count
                else:
                    # Legacy mode - access the merged file
                    _ = self.host_fasta
                    stats['fasta']['hosts'] = self._host_count if self._host_count else len(self.host_fasta.keys())
            except Exception as e:
                logging.warning(f"Host data configured but tables not found. Run host genome workflow first. Error: {e}")
                self._has_host_data = False  # Disable host support if tables don't exist

        try:
            source_breakdown = self.conn.execute(
                """
                SELECT source_type, Source_DB, COUNT(*) AS count
                FROM fact_phages
                GROUP BY source_type, Source_DB
                ORDER BY source_type, count DESC
                """
            ).fetchdf()
            stats['database']['source_breakdown'] = source_breakdown.to_dict(orient='records')
        except Exception as e:
            logging.debug(f"Could not compute source breakdown stats: {e}")
        
        logging.info(f"📊 Database Stats:")
        logging.info(f"   Phages: {stats['database']['phages']:,}")
        logging.info(f"   Proteins: {stats['database']['proteins']:,}")
        if 'hosts' in stats['database']:
            logging.info(f"   Hosts: {stats['database']['hosts']:,}")
            if 'phage_host_associations' in stats['database']:
                logging.info(f"   Phage-Host Associations: {stats['database']['phage_host_associations']:,}")
        
        logging.info(f"📊 FASTA Stats:")
        logging.info(f"   Phages: {stats['fasta']['phages']:,}")
        logging.info(f"   Proteins: {stats['fasta']['proteins']:,}")
        if 'hosts' in stats['fasta']:
            logging.info(f"   Hosts: {stats['fasta']['hosts']:,}")
        
        return stats

    def get_phages(
        self,
        source_type: Optional[str] = None,
        source_db: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """Return phage metadata with optional provenance filters."""
        filters = []
        params = []
        if source_type:
            filters.append("source_type = ?")
            params.append(source_type)
        if source_db:
            filters.append("Source_DB = ?")
            params.append(source_db)

        query = "SELECT * FROM fact_phages"
        if filters:
            query += " WHERE " + " AND ".join(filters)
        query += " ORDER BY Phage_ID"
        if limit:
            query += f" LIMIT {int(limit)}"

        return self.conn.execute(query, params).fetchdf()

    def get_interactions(
        self,
        source_type: Optional[str] = None,
        source_db: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """Return phage-host interactions with optional provenance filters."""
        if not self._has_host_data:
            raise ValueError("Host data not available - run host genome download workflow first")

        filters = []
        params = []
        if source_type:
            filters.append("p.source_type = ?")
            params.append(source_type)
        if source_db:
            filters.append("p.Source_DB = ?")
            params.append(source_db)

        query = """
        SELECT DISTINCT
            pha.Phage_ID,
            pha.Host_ID,
            p.Source_DB,
            CASE
                WHEN LOWER(TRIM(COALESCE(p.source_type, ''))) = 'private' THEN 'private'
                ELSE 'public'
            END as source_type
        FROM phage_host_associations pha
        JOIN fact_phages p ON p.Phage_ID = pha.Phage_ID
        """
        if filters:
            query += " WHERE " + " AND ".join(filters)
        query += " ORDER BY pha.Phage_ID, pha.Host_ID"
        if limit:
            query += f" LIMIT {int(limit)}"

        return self.conn.execute(query, params).fetchdf()

    def get_pair_ids_only(
        self,
        shuffle: bool = False,
        seed: Optional[int] = None,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """Return distinct (Phage_ID, Host_ID) pairs without loading sequences.

        Lightweight alternative to :meth:`get_phage_host_pairs` when only the
        pair identifiers are needed (e.g. for counting, sampling, or building
        negative-example sets).

        Args:
            shuffle: If ``True``, order rows deterministically using the MD5
                hash of the concatenated pair IDs.  This produces a stable,
                pseudo-random ordering without requiring Python-side
                shuffling, which is useful for reproducible sub-sampling of
                large pair lists.
            seed: Reserved for future use.  The MD5 ordering is inherently
                deterministic so *seed* currently has no effect.
            limit: Optional maximum number of pairs to return.

        Returns:
            DataFrame with columns ``Phage_ID`` and ``Host_ID``.

        Example:
            >>> # All pairs in deterministic (alphabetical) order
            >>> pairs = retriever.get_pair_ids_only()
            >>> len(pairs)
            1240905

            >>> # Deterministic pseudo-random order for reproducible sampling
            >>> pairs = retriever.get_pair_ids_only(shuffle=True)
            >>> pairs.head()

            >>> # First 1000 pairs in shuffled order
            >>> sample = retriever.get_pair_ids_only(shuffle=True, limit=1000)
        """
        query = "SELECT DISTINCT Phage_ID, Host_ID FROM phage_host_associations"
        if shuffle:
            query += " ORDER BY MD5(Phage_ID || Host_ID)"
        if limit:
            query += f" LIMIT {int(limit)}"
        return self.conn.execute(query).fetchdf()

    def get_host_sequences(self, query: str, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get host sequences based on SQL query
        
        Args:
            query: SQL query that returns Host_ID column
            limit: Optional limit on number of sequences
        
        Returns:
            DataFrame with columns: Host_ID, Species_Name, Sequence, Length, GC_Content
        
        Example:
            query = "SELECT Host_ID FROM dim_hosts WHERE Species_Name LIKE '%Escherichia%'"
            df = retriever.get_host_sequences(query)
        """
        if not self._has_host_data:
            raise ValueError(
                "Host data not available. Please run the host genome download workflow first:\n"
                "  snakemake --use-conda --cores 1 all_hosts\n"
                "Or check that host_fasta_path or host_mapping_path was provided when creating SequenceRetriever."
            )
        
        # In legacy mode, ensure host FASTA is loaded
        if not self._use_host_mapping:
            _ = self.host_fasta
        
        logging.info(f"🔍 Executing query: {query[:100]}...")
        
        if limit:
            query = f"{query} LIMIT {limit}"
        
        result = self.conn.execute(query).fetchdf()
        
        if 'Host_ID' not in result.columns:
            raise ValueError("Query must return 'Host_ID' column")
        
        host_ids = result['Host_ID'].tolist()
        logging.info(f"📊 Retrieved {len(host_ids):,} Host IDs from query")
        
        return self._fetch_host_sequences(host_ids)
    
    def get_host_by_phage(self, phage_id: str) -> pd.DataFrame:
        """
        Get host genome(s) for a given phage
        
        Args:
            phage_id: Phage ID
        
        Returns:
            DataFrame with host sequences associated with the phage
        
        Example:
            df = retriever.get_host_by_phage("NC_000866")
        """
        if not self._has_host_data:
            raise ValueError("Host data not available - run host genome download workflow first")
        
        query = f"""
        SELECT h.Host_ID
        FROM phage_host_associations pha
        JOIN dim_hosts h ON pha.Host_ID = h.Host_ID
        WHERE pha.Phage_ID = '{phage_id}'
        """
        
        return self.get_host_sequences(query)
    
    def get_phage_host_pairs(
        self,
        where_clause: str = None,
        limit: Optional[int] = None,
        host_contig_mode: str = "concat",
        phage_contig_mode: str = "first",
    ) -> pd.DataFrame:
        """
        Get phage-host interaction pairs with sequences and metadata.

        Args:
            where_clause: Optional SQL WHERE clause to filter pairs (without
                the ``WHERE`` keyword).
            limit: Optional limit on number of pairs.
            host_contig_mode: How to handle multi-contig host FASTA files.

                * ``"first"`` – return only the first/largest
                  contig (original behaviour; fully backward-compatible).
                * ``"concat"`` *(default)* – concatenate all contigs into one string
                  (sorted by length desc).  Use this to include the full
                  host genome even when it is split across scaffolds.

            phage_contig_mode: How to handle phage sequences.  Currently
                phage FASTA records are always single-contig, so both
                ``"first"`` and ``"concat"`` are equivalent.  Provided for
                forward-compatibility.

        Returns:
            DataFrame with columns: Phage_ID, Host_ID, Phage_Source,
            Phage_Length, Phage_GC, Phage_Taxonomy, Phage_Completeness,
            Phage_Lifestyle, Phage_Cluster, Phage_Subcluster, Species_Name,
            Host_Assembly_Level, Host_Length, Host_GC, Host_RefSeq_Category,
            Phage_Sequence, Host_Sequence.

        Example:
            # Get all pairs (default – full host genome, concatenated contigs)
            pairs = retriever.get_phage_host_pairs()

            # Full host genome even for fragmented assemblies
            pairs = retriever.get_phage_host_pairs(host_contig_mode="concat")

            # Full host genome with 100-N gaps between scaffolds
            # (use get_phage_host_pairs_iterator for gap control)

            # Get pairs for specific lifestyle
            pairs = retriever.get_phage_host_pairs("p.Lifestyle = 'Lytic'", limit=1000)

            # Get pairs with complete host genomes
            pairs = retriever.get_phage_host_pairs("h.Assembly_Level = 'Complete Genome'")
        """
        if not self._has_host_data:
            raise ValueError("Host data not available - run host genome download workflow first")

        # Validate contig_mode arguments early
        for mode_name, mode_val in (
            ("host_contig_mode", host_contig_mode),
            ("phage_contig_mode", phage_contig_mode),
        ):
            if mode_val not in ("first", "concat"):
                raise ValueError(
                    f"{mode_name} must be 'first' or 'concat', got '{mode_val}'."
                )

        # Ensure public phage FASTA is loaded (only strictly needed for public phages,
        # but we load it eagerly since most datasets have at least some public phages).
        _ = self.phage_fasta
        # For legacy host mode, ensure host FASTA is loaded
        if not self._use_host_mapping:
            _ = self.host_fasta

        # Build query
        query = """
        SELECT DISTINCT
            pha.Phage_ID,
            pha.Host_ID,
            p.Source_DB as Phage_Source,
            CASE
                WHEN LOWER(TRIM(COALESCE(p.source_type, ''))) = 'private' THEN 'private'
                ELSE 'public'
            END as Phage_Source_Type,
            p.Length as Phage_Length,
            p.GC_content as Phage_GC,
            p.Taxonomy as Phage_Taxonomy,
            p.Completeness as Phage_Completeness,
            p.Lifestyle as Phage_Lifestyle,
            p.Cluster as Phage_Cluster,
            p.Subcluster as Phage_Subcluster,
            h.Species_Name,
            h.Assembly_Level as Host_Assembly_Level,
            h.Genome_Length as Host_Length,
            h.GC_Content as Host_GC,
            h.RefSeq_Category as Host_RefSeq_Category
        FROM phage_host_associations pha
        JOIN fact_phages p ON pha.Phage_ID = p.Phage_ID
        JOIN dim_hosts h ON pha.Host_ID = h.Host_ID
        """

        # Parse where_clause to separate WHERE conditions from LIMIT/OFFSET
        where_conditions, limit_offset = parse_where_clause(where_clause)

        if where_conditions:
            query += f" WHERE {where_conditions}"

        # If limit parameter is provided, it takes precedence over any LIMIT in where_clause
        if limit:
            query += f" LIMIT {limit}"
        elif limit_offset:
            query += f" {limit_offset}"

        logging.info(f"🔍 Querying phage-host pairs...")
        result = self.conn.execute(query).fetchdf()

        logging.info(f"📊 Found {len(result):,} phage-host pairs")

        # Fetch sequences
        phage_ids = result['Phage_ID'].tolist()
        host_ids = result['Host_ID'].tolist()

        logging.info(
            f"📥 Fetching sequences for {len(phage_ids):,} phages and "
            f"{len(set(host_ids)):,} unique hosts"
        )

        phage_seqs = {}
        host_seqs = {}

        # Build per-phage source lookup so private phages are routed to their
        # per-source FASTA rather than all_phages.fasta.
        phage_source_db = dict(zip(result['Phage_ID'], result['Phage_Source']))
        phage_source_type = dict(zip(result['Phage_ID'], result['Phage_Source_Type']))

        # Fetch phage sequences (routes public → all_phages.fasta, private → per-source FASTA)
        for phage_id in phage_ids:
            seq = self._get_phage_sequence(
                phage_id,
                source_db=phage_source_db.get(phage_id),
                source_type=phage_source_type.get(phage_id),
            )
            phage_seqs[phage_id] = seq

        # Fetch host sequences (unique only)
        for host_id in set(host_ids):
            try:
                seq = self.get_host_sequence(host_id, contig_mode=host_contig_mode)
                host_seqs[host_id] = seq
            except KeyError:
                host_seqs[host_id] = None

        # Add sequences to result
        result['Phage_Sequence'] = result['Phage_ID'].map(phage_seqs)
        result['Host_Sequence'] = result['Host_ID'].map(host_seqs)

        missing_phage_ids = sorted(
            result.loc[result['Phage_Sequence'].isna(), 'Phage_ID'].drop_duplicates().tolist()
        )
        missing_host_ids = sorted(
            result.loc[result['Host_Sequence'].isna(), 'Host_ID'].drop_duplicates().tolist()
        )

        # Filter out rows with missing sequences
        before_count = len(result)
        result = result.dropna(subset=['Phage_Sequence', 'Host_Sequence'])
        after_count = len(result)

        if before_count > after_count:
            logging.warning(
                "⚠️  Removed %d pairs with missing sequences "
                "(%d phages missing, %d hosts missing)",
                before_count - after_count,
                len(missing_phage_ids),
                len(missing_host_ids),
            )
            if missing_phage_ids:
                logging.warning("   Missing phage IDs (sample): %s", ", ".join(missing_phage_ids[:5]))
            if missing_host_ids:
                logging.warning("   Missing host IDs (sample): %s", ", ".join(missing_host_ids[:5]))

        logging.info(f"✅ Retrieved {len(result):,} complete phage-host pairs with sequences")

        return result

    def query_phage_host_pairs(
        self,
        phage_filters: Optional[Dict[str, str]] = None,
        host_filters: Optional[Dict[str, str]] = None,
        limit: Optional[int] = None,
        host_contig_mode: str = "concat",
        phage_contig_mode: str = "first",
    ) -> pd.DataFrame:
        """
        Query phage-host pairs using structured filter dictionaries.

        A convenience wrapper around :meth:`get_phage_host_pairs` that accepts
        column-level filter dicts instead of a raw SQL WHERE clause.

        Args:
            phage_filters: Mapping of phage column names to filter values.
                Keys are columns in the ``fact_phages`` table (e.g.
                ``"Lifestyle"``, ``"Source_DB"``).  Values that contain ``%``
                are matched with SQL ``LIKE``; all other values use
                case-insensitive equality (``LOWER(col) = LOWER(val)``).
            host_filters: Mapping of host column names to filter values.
                Keys are columns in the ``dim_hosts`` table.
                ``"Organism_Name"`` is accepted as an alias for
                ``"Species_Name"``.  The same ``%`` / equality rule applies.
            limit: Maximum number of pairs to return.
            host_contig_mode: See :meth:`get_phage_host_pairs`.
            phage_contig_mode: See :meth:`get_phage_host_pairs`.

        Returns:
            DataFrame — same schema as :meth:`get_phage_host_pairs`.

        Example:
            # Lytic phages infecting Escherichia hosts
            pairs = retriever.query_phage_host_pairs(
                phage_filters={'Lifestyle': 'lytic'},
                host_filters={'Organism_Name': '%Escherichia%'},
                limit=20,
            )

            # Complete-genome hosts only
            pairs = retriever.query_phage_host_pairs(
                host_filters={'Assembly_Level': 'Complete Genome'},
            )
        """
        # Allowed column names for each table (guards against column-name injection)
        _ALLOWED_PHAGE_COLS: frozenset = frozenset({
            "Phage_ID", "Source_DB", "Length", "GC_content", "Taxonomy",
            "Completeness", "Host", "Lifestyle", "Cluster", "Subcluster",
            "source_type",
        })
        _ALLOWED_HOST_COLS: frozenset = frozenset({
            "Host_ID", "Species_Name", "Strain_Name", "Assembly_Accession",
            "Assembly_Name", "Assembly_Level", "Genome_Length", "GC_Content",
            "RefSeq_Category", "Download_Date", "Source",
        })
        # Map intuitive host column aliases to the actual dim_hosts column names
        HOST_COLUMN_ALIASES: Dict[str, str] = {
            "Organism_Name": "Species_Name",
        }

        conditions: List[str] = []

        for col, val in (phage_filters or {}).items():
            if col not in _ALLOWED_PHAGE_COLS:
                raise ValueError(
                    f"Unknown phage filter column '{col}'. "
                    f"Allowed columns: {sorted(_ALLOWED_PHAGE_COLS)}"
                )
            col_expr = f"p.{col}"
            val_str = str(val).replace("'", "''")
            if "%" in val_str:
                conditions.append(f"{col_expr} LIKE '{val_str}'")
            else:
                conditions.append(f"LOWER({col_expr}) = LOWER('{val_str}')")

        for col, val in (host_filters or {}).items():
            real_col = HOST_COLUMN_ALIASES.get(col, col)
            if real_col not in _ALLOWED_HOST_COLS:
                raise ValueError(
                    f"Unknown host filter column '{col}'. "
                    f"Allowed columns: {sorted(_ALLOWED_HOST_COLS | set(HOST_COLUMN_ALIASES))}"
                )
            col_expr = f"h.{real_col}"
            val_str = str(val).replace("'", "''")
            if "%" in val_str:
                conditions.append(f"{col_expr} LIKE '{val_str}'")
            else:
                conditions.append(f"LOWER({col_expr}) = LOWER('{val_str}')")

        where_clause = " AND ".join(conditions) if conditions else None

        return self.get_phage_host_pairs(
            where_clause=where_clause,
            limit=limit,
            host_contig_mode=host_contig_mode,
            phage_contig_mode=phage_contig_mode,
        )

    def get_phage_sequence(self, phage_id: str) -> Optional[str]:
        """
        Get the DNA sequence for a single phage.

        A public convenience wrapper around the internal
        :meth:`_get_phage_sequence` method for use in notebooks and scripts
        where the source database is not known in advance.  The method first
        tries the public ``all_phages.fasta`` and falls back to private
        per-source FASTA files when private phage mapping is configured.

        Args:
            phage_id: Phage identifier (first token of the FASTA header).

        Returns:
            Sequence string, or ``None`` if the phage is not found in any
            indexed FASTA.

        Example:
            seq = retriever.get_phage_sequence('PHAGE_Esche_phiX174_NC_001422')
            if seq:
                print(f'Length: {len(seq):,} bp')
        """
        return self._get_phage_sequence(phage_id)

    def _fetch_host_sequences(self, host_ids: list) -> pd.DataFrame:
        """
        Fetch host sequences for given Host IDs
        
        Args:
            host_ids: List of Host IDs to retrieve sequences for
        
        Returns:
            DataFrame with Host_ID and Sequence columns
        """
        if not host_ids:
            logging.warning("No Host IDs provided")
            return pd.DataFrame(columns=['Host_ID', 'Sequence', 'Length', 'GC_Content'])
        
        logging.info(f"🔍 Fetching sequences for {len(host_ids):,} hosts")
        
        sequences = []
        missing_ids = []
        
        for host_id in host_ids:
            try:
                seq_str = self.get_host_sequence(host_id)
                sequences.append({
                    'Host_ID': host_id,
                    'Sequence': seq_str,
                    'Length': len(seq_str),
                    'GC_Content': round(gc_fraction(seq_str) * 100, 2) if len(seq_str) > 0 else 0.0
                })
            except KeyError:
                missing_ids.append(host_id)
                logging.warning(f"⚠️  Host ID '{host_id}' not found in FASTA")
        
        if missing_ids:
            logging.warning(f"⚠️  {len(missing_ids):,} host IDs not found in FASTA file")
        
        df = pd.DataFrame(sequences)
        logging.info(f"✅ Retrieved {len(df):,} sequences")
        
        return df
    
    def get_phage_metadata(self, where_clause: str = None, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get phage metadata from the database
        
        Args:
            where_clause: Optional SQL WHERE clause to filter phages (without WHERE keyword)
            limit: Optional limit on number of phages
        
        Returns:
            DataFrame with phage metadata including: Phage_ID, Source_DB, Length, GC_content,
            Taxonomy, Completeness, Host, Lifestyle, Cluster, Subcluster
        
        Example:
            # Get all phages metadata
            metadata = retriever.get_phage_metadata()
            
            # Get phages from specific source
            metadata = retriever.get_phage_metadata("Source_DB = 'PhagesDB'", limit=1000)
            
            # Get lytic phages
            metadata = retriever.get_phage_metadata("Lifestyle = 'Lytic'")
        """
        query = """
        SELECT 
            Phage_ID,
            Source_DB,
            Length,
            GC_content,
            Taxonomy,
            Completeness,
            Host,
            Lifestyle,
            Cluster,
            Subcluster
        FROM fact_phages
        """
        
        # Parse where_clause to separate WHERE conditions from LIMIT/OFFSET
        where_conditions, limit_offset = parse_where_clause(where_clause)
        
        if where_conditions:
            query += f" WHERE {where_conditions}"
        
        # If limit parameter is provided, it takes precedence over any LIMIT in where_clause
        if limit:
            query += f" LIMIT {limit}"
        elif limit_offset:
            query += f" {limit_offset}"
        
        logging.info(f"🔍 Querying phage metadata...")
        result = self.conn.execute(query).fetchdf()
        
        logging.info(f"✅ Retrieved metadata for {len(result):,} phages")
        
        return result
    
    def get_host_metadata(self, where_clause: str = None, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get host metadata from the database
        
        Args:
            where_clause: Optional SQL WHERE clause to filter hosts (without WHERE keyword)
            limit: Optional limit on number of hosts
        
        Returns:
            DataFrame with host metadata including: Host_ID, Species_Name, Strain_Name,
            Assembly_Accession, Assembly_Name, Assembly_Level, Genome_Length, GC_Content,
            RefSeq_Category, Download_Date, Source
        
        Example:
            # Get all hosts metadata
            metadata = retriever.get_host_metadata()
            
            # Get hosts of specific species
            metadata = retriever.get_host_metadata("Species_Name LIKE '%Escherichia%'")
            
            # Get complete genomes only
            metadata = retriever.get_host_metadata("Assembly_Level = 'Complete Genome'")
        """
        if not self._has_host_data:
            raise ValueError("Host data not available - run host genome download workflow first")
        
        query = """
        SELECT 
            Host_ID,
            Species_Name,
            Strain_Name,
            Assembly_Accession,
            Assembly_Name,
            Assembly_Level,
            Genome_Length,
            GC_Content,
            RefSeq_Category,
            Download_Date,
            Source
        FROM dim_hosts
        """
        
        # Parse where_clause to separate WHERE conditions from LIMIT/OFFSET
        where_conditions, limit_offset = parse_where_clause(where_clause)
        
        if where_conditions:
            query += f" WHERE {where_conditions}"
        
        # If limit parameter is provided, it takes precedence over any LIMIT in where_clause
        if limit:
            query += f" LIMIT {limit}"
        elif limit_offset:
            query += f" {limit_offset}"
        
        logging.info(f"🔍 Querying host metadata...")
        result = self.conn.execute(query).fetchdf()
        
        logging.info(f"✅ Retrieved metadata for {len(result):,} hosts")
        
        return result
    
    def get_phage_host_metadata(self, where_clause: str = None, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get combined phage-host metadata for interaction pairs
        
        Args:
            where_clause: Optional SQL WHERE clause to filter pairs (without WHERE keyword)
            limit: Optional limit on number of pairs
        
        Returns:
            DataFrame with combined phage and host metadata
        
        Example:
            # Get all pairs metadata
            metadata = retriever.get_phage_host_metadata()
            
            # Get pairs from specific phage source
            metadata = retriever.get_phage_host_metadata("p.Source_DB = 'PhagesDB'")
            
            # Get pairs with lytic phages only
            metadata = retriever.get_phage_host_metadata("p.Lifestyle = 'Lytic'", limit=1000)
        """
        if not self._has_host_data:
            raise ValueError("Host data not available - run host genome download workflow first")
        
        query = """
        SELECT DISTINCT
            pha.Phage_ID,
            pha.Host_ID,
            p.Source_DB as Phage_Source,
            CASE
                WHEN LOWER(TRIM(COALESCE(p.source_type, ''))) = 'private' THEN 'private'
                ELSE 'public'
            END as Phage_Source_Type,
            p.Length as Phage_Length,
            p.GC_content as Phage_GC,
            p.Taxonomy as Phage_Taxonomy,
            p.Completeness as Phage_Completeness,
            p.Lifestyle as Phage_Lifestyle,
            p.Cluster as Phage_Cluster,
            p.Subcluster as Phage_Subcluster,
            h.Species_Name as Host_Species,
            h.Strain_Name as Host_Strain,
            h.Assembly_Accession as Host_Assembly,
            h.Assembly_Level as Host_Assembly_Level,
            h.Genome_Length as Host_Length,
            h.GC_Content as Host_GC,
            h.RefSeq_Category as Host_RefSeq_Category,
            h.Source as Host_Source
        FROM phage_host_associations pha
        JOIN fact_phages p ON pha.Phage_ID = p.Phage_ID
        JOIN dim_hosts h ON pha.Host_ID = h.Host_ID
        """
        
        # Parse where_clause to separate WHERE conditions from LIMIT/OFFSET
        where_conditions, limit_offset = parse_where_clause(where_clause)
        
        if where_conditions:
            query += f" WHERE {where_conditions}"
        
        # If limit parameter is provided, it takes precedence over any LIMIT in where_clause
        if limit:
            query += f" LIMIT {limit}"
        elif limit_offset:
            query += f" {limit_offset}"
        
        logging.info(f"🔍 Querying phage-host metadata...")
        result = self.conn.execute(query).fetchdf()
        
        logging.info(f"✅ Retrieved metadata for {len(result):,} phage-host pairs")
        
        return result
    
    def get_protein_metadata(self, where_clause: str = None, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get protein metadata from dim_proteins
        
        Args:
            where_clause: Optional SQL WHERE clause to filter proteins (without WHERE keyword)
            limit: Optional limit on number of proteins
        
        Returns:
            DataFrame with protein metadata
        
        Example:
            # Get all proteins metadata
            metadata = retriever.get_protein_metadata()
            
            # Get proteins from specific source
            metadata = retriever.get_protein_metadata("Source_DB = 'PhagesDB'", limit=1000)
        """
        query = """
        SELECT 
            Protein_ID,
            Phage_ID,
            Source_DB,
            Protein_source,
            Function_prediction_source,
            Start,
            Stop,
            Strand,
            Product,
            Protein_classification
        FROM dim_proteins
        """
        
        # Parse where_clause to separate WHERE conditions from LIMIT/OFFSET
        where_conditions, limit_offset = parse_where_clause(where_clause)
        
        if where_conditions:
            query += f" WHERE {where_conditions}"
        
        # If limit parameter is provided, it takes precedence over any LIMIT in where_clause
        if limit:
            query += f" LIMIT {limit}"
        elif limit_offset:
            query += f" {limit_offset}"
        
        logging.info(f"🔍 Querying protein metadata...")
        result = self.conn.execute(query).fetchdf()
        
        logging.info(f"✅ Retrieved metadata for {len(result):,} proteins")
        
        return result
    
    def help(self):
        """Print help information"""
        help_text = """
        SequenceRetriever Help:

        Core sequence retrieval methods:
            - get_phage_sequences(query: str, limit: Optional[int] = None) -> pd.DataFrame
            - get_protein_sequences(query: str, limit: Optional[int] = None) -> pd.DataFrame
            - get_host_sequences(query: str, limit: Optional[int] = None) -> pd.DataFrame
            - get_host_sequence(host_id: str, contig_mode: str = "first") -> str
                  contig_mode="first"   : return only the first/largest contig (default)
                  contig_mode="concat"  : concatenate all contigs into one string

        Full-genome retrieval (multi-contig support):
            - get_host_genome(host_id, mode="concat", gap=0, order="length_desc")
                  Retrieve complete host genome even when split across scaffolds.
                  mode="concat"  : all contigs joined into one string (default)
                  mode="first"   : only the first/largest contig
                  mode="list"    : list of per-contig strings
                  mode="dict"    : {header: sequence} mapping
                  gap=N          : insert N "N" characters between contigs
                  order="length_desc" : sort by length desc (default, deterministic)
                  order="file"        : preserve FASTA file order
            - get_host_genome_stats(host_id, order="length_desc") -> dict
                  Returns {"contig_count", "lengths", "total_length"}
            - get_phage_genome(phage_id, mode="concat", gap=0, order="length_desc")
                  Same interface as get_host_genome; for phages (usually single-contig).

        Pair retrieval methods:
            - get_phage_host_pairs(where_clause=None, limit=None,
                                   host_contig_mode="concat", phage_contig_mode="first")
                   -> pd.DataFrame  with Phage_Sequence, Host_Sequence columns
                   Default host mode returns full fragmented host genomes.
            - get_phage_host_pairs_iterator(where_clause=None, batch_size=1000,
                                            host_contig_mode="first", phage_contig_mode="first")
                  -> iterator of DataFrame batches (memory-efficient)

        Metadata methods:
            - get_sequences_by_ids(phage_ids, protein_ids) -> Dict
            - get_protein_sequences_by_phage(phage_id: str) -> pd.DataFrame
            - get_phage_metadata(where_clause=None, limit=None) -> pd.DataFrame
            - get_host_metadata(where_clause=None, limit=None) -> pd.DataFrame
            - get_phage_host_metadata(where_clause=None, limit=None) -> pd.DataFrame
            - export_fasta(df, output_path, id_col="Phage_ID")
            - get_stats() -> Dict
            - close()

        Usage Examples:
            # Connect
            retriever = SequenceRetriever(db_path, phage_fasta_path, protein_fasta_path,
                                          host_mapping_path=host_mapping_path)

            # Standard sequence retrieval (unchanged)
            phage_df = retriever.get_phage_sequences(
                "SELECT Phage_ID FROM fact_phages WHERE Length > 50000", limit=100)

            # Get full host genome (multi-contig safe)
            full_genome = retriever.get_host_genome("GCF_000005845")
            print(len(full_genome), "bp total")

            # Inspect contig fragmentation
            stats = retriever.get_host_genome_stats("GCF_000005845")
            print(stats["contig_count"], "contigs,", stats["total_length"], "bp total")

            # Phage-host pairs with concatenated host genomes
            pairs = retriever.get_phage_host_pairs(
                "p.Lifestyle = 'Lytic'", limit=100, host_contig_mode="concat")

            # Batch iterator with full host genomes
            for batch_df in retriever.get_phage_host_pairs_iterator(
                    host_contig_mode="concat", batch_size=500):
                print(f"Batch: {len(batch_df)} pairs, "
                      f"host genome sizes: {batch_df['Host_Sequence'].str.len().describe()}")

            # Export, stats, close
            retriever.export_fasta(phage_df, "output_phages.fasta", id_col="Phage_ID")
            retriever.get_stats()
            retriever.close()
        """
        print(help_text)
    
    def create_streaming_dataset(
        self,
        where_clause: Optional[str] = None,
        batch_size: int = 1000,
        transform: Optional[object] = None,
        missing_hosts_csv: Optional[str] = None
    ):
        """
        Create a PhageHostStreamingDataset for memory-efficient iteration.
        
        This factory method creates a streaming dataset that fetches data in batches
        from DuckDB and loads sequences on-demand. Ideal for large datasets.
        
        Args:
            where_clause: Optional SQL WHERE clause to filter pairs (without WHERE keyword)
            batch_size: Number of records to fetch per database query (default: 1000)
            transform: Optional transform function to apply to each sample
            missing_hosts_csv: Optional path to save CSV of phages with missing hosts
                              (e.g., "/data/intermediate/missing_hosts.csv")
            
        Returns:
            PhageHostStreamingDataset instance
            
        Example:
            >>> dataset = retriever.create_streaming_dataset(
            ...     where_clause="Confidence > 0.8",
            ...     batch_size=1000,
            ...     missing_hosts_csv="/data/intermediate/missing_hosts.csv"
            ... )
            >>> from torch.utils.data import DataLoader
            >>> dataloader = DataLoader(dataset, batch_size=32)
            >>> for batch in dataloader:
            ...     # Process batch
            ...     pass
        """
        from .streaming_dataset import PhageHostStreamingDataset
        
        # Get database path from connection
        # PRAGMA database_list returns (seq, name, file) - index 2 is the file path
        db_path = str(self.conn.execute("PRAGMA database_list").fetchone()[2])
        
        return PhageHostStreamingDataset(
            db_path=db_path,
            phage_fasta_path=self._phage_fasta_path,
            host_fasta_path=self._host_fasta_path,
            host_mapping_path=self._host_mapping_path,
            where_clause=where_clause,
            batch_size=batch_size,
            transform=transform,
            missing_hosts_csv=missing_hosts_csv
        )
    
    def create_indexed_dataset(
        self,
        where_clause: Optional[str] = None,
        transform: Optional[object] = None,
        missing_hosts_csv: Optional[str] = None
    ):
        """
        Create a PhageHostIndexedDataset for random access with caching.
        
        This factory method creates an indexed dataset that caches metadata in memory
        and provides random access. Suitable for medium-sized datasets.
        
        Args:
            where_clause: Optional SQL WHERE clause to filter pairs (without WHERE keyword)
            transform: Optional transform function to apply to each sample
            missing_hosts_csv: Optional path to save CSV of phages with missing hosts
                              (e.g., "/data/intermediate/missing_hosts.csv")
            
        Returns:
            PhageHostIndexedDataset instance
            
        Example:
            >>> dataset = retriever.create_indexed_dataset(
            ...     where_clause="Confidence > 0.8",
            ...     missing_hosts_csv="/data/intermediate/missing_hosts.csv"
            ... )
            >>> from torch.utils.data import DataLoader
            >>> # Supports shuffling and multi-worker loading
            >>> dataloader = DataLoader(dataset, batch_size=32, shuffle=True, num_workers=4)
            >>> for batch in dataloader:
            ...     # Process batch
            ...     pass
        """
        from .streaming_dataset import PhageHostIndexedDataset
        
        # Get database path from connection
        # PRAGMA database_list returns (seq, name, file) - index 2 is the file path
        db_path = str(self.conn.execute("PRAGMA database_list").fetchone()[2])
        
        return PhageHostIndexedDataset(
            db_path=db_path,
            phage_fasta_path=self._phage_fasta_path,
            host_fasta_path=self._host_fasta_path,
            host_mapping_path=self._host_mapping_path,
            where_clause=where_clause,
            transform=transform,
            missing_hosts_csv=missing_hosts_csv
        )
    
    def get_phage_host_pairs_iterator(
        self,
        where_clause: Optional[str] = None,
        batch_size: int = 1000,
        host_contig_mode: str = "first",
        phage_contig_mode: str = "first",
    ):
        """Get an iterator that yields batches of phage-host pairs as DataFrames.

        This provides a simple memory-efficient way to process large datasets
        in batches without loading everything into memory at once.
        Alternative for non-PyTorch workflows.

        Args:
            where_clause: Optional SQL WHERE clause to filter pairs (without
                the ``WHERE`` keyword).
            batch_size: Number of pairs to fetch per batch (default: 1000).
            host_contig_mode: How to handle multi-contig host FASTA files.

                * ``"first"`` *(default)* – return only the first/largest
                  contig (original behaviour; fully backward-compatible).
                * ``"concat"`` – concatenate all contigs into one string
                  (sorted by length desc).  Use this to include the full
                  host genome even when it is split across scaffolds.

            phage_contig_mode: How to handle phage sequences.  Currently
                phage FASTA records are always single-contig, so both
                ``"first"`` and ``"concat"`` are equivalent.  Provided for
                forward-compatibility.

        Yields:
            DataFrame batches containing phage-host pairs with sequences and
            metadata.  The ``Host_Sequence`` column will contain the
            concatenated genome when *host_contig_mode* is ``"concat"``.

        Example:
            >>> # Default – single contig per host (original behaviour)
            >>> for batch_df in retriever.get_phage_host_pairs_iterator(
            ...     where_clause="Confidence > 0.8",
            ...     batch_size=1000,
            ... ):
            ...     print(f"Processing {len(batch_df)} pairs")

            >>> # Full concatenated host genome
            >>> for batch_df in retriever.get_phage_host_pairs_iterator(
            ...     host_contig_mode="concat",
            ...     batch_size=500,
            ... ):
            ...     print(batch_df[["Phage_ID", "Host_ID", "Host_Sequence"]].head())
        """
        if not self._has_host_data:
            raise ValueError("Host data not available - run host genome download workflow first")

        # Validate contig_mode arguments early
        for mode_name, mode_val in (
            ("host_contig_mode", host_contig_mode),
            ("phage_contig_mode", phage_contig_mode),
        ):
            if mode_val not in ("first", "concat"):
                raise ValueError(
                    f"{mode_name} must be 'first' or 'concat', got '{mode_val}'."
                )

        # Ensure FASTA files are loaded
        _ = self.phage_fasta
        if not self._use_host_mapping:
            _ = self.host_fasta

        # Build query
        query = """
        SELECT DISTINCT
            pha.Phage_ID,
            pha.Host_ID,
            p.Source_DB as Phage_Source,
            COALESCE(NULLIF(NULLIF(LOWER(TRIM(p.source_type)), 'nan'), ''), 'public') as Phage_Source_Type,
            p.Length as Phage_Length,
            p.GC_content as Phage_GC,
            p.Taxonomy as Phage_Taxonomy,
            p.Completeness as Phage_Completeness,
            p.Lifestyle as Phage_Lifestyle,
            p.Cluster as Phage_Cluster,
            p.Subcluster as Phage_Subcluster,
            h.Species_Name,
            h.Assembly_Level as Host_Assembly_Level,
            h.Genome_Length as Host_Length,
            h.GC_Content as Host_GC,
            h.RefSeq_Category as Host_RefSeq_Category
        FROM phage_host_associations pha
        JOIN fact_phages p ON pha.Phage_ID = p.Phage_ID
        JOIN dim_hosts h ON pha.Host_ID = h.Host_ID
        """

        # Parse where_clause to separate WHERE conditions from LIMIT/OFFSET
        where_conditions, limit_offset = parse_where_clause(where_clause)

        if where_conditions:
            query += f" WHERE {where_conditions}"

        if limit_offset:
            query += f" {limit_offset}"

        logging.info(f"🔍 Starting batch iteration with batch_size={batch_size}")

        # Execute query and fetch in batches
        cursor = self.conn.execute(query)
        batch_num = 0

        while True:
            # Fetch a batch
            batch_df = cursor.fetch_df_chunk(batch_size)
            if batch_df is None or len(batch_df) == 0:
                break

            batch_num += 1
            logging.info(f"📦 Processing batch {batch_num} ({len(batch_df)} pairs)")

            # Fetch sequences for this batch
            phage_seqs = {}
            host_seqs = {}

            phage_source_db = dict(zip(batch_df['Phage_ID'], batch_df['Phage_Source']))
            phage_source_type = dict(zip(batch_df['Phage_ID'], batch_df['Phage_Source_Type']))

            for phage_id in batch_df['Phage_ID'].unique():
                phage_seqs[phage_id] = self._get_phage_sequence(
                    phage_id,
                    source_db=phage_source_db.get(phage_id),
                    source_type=phage_source_type.get(phage_id),
                )

            for host_id in batch_df['Host_ID'].unique():
                host_seqs[host_id] = self._get_sequence_safe(
                    host_id, 'host', host_contig_mode=host_contig_mode
                )

            # Add sequences to batch
            batch_df['Phage_Sequence'] = batch_df['Phage_ID'].map(phage_seqs)
            batch_df['Host_Sequence'] = batch_df['Host_ID'].map(host_seqs)

            missing_phage_ids = sorted(
                batch_df.loc[batch_df['Phage_Sequence'].isna(), 'Phage_ID'].drop_duplicates().tolist()
            )
            missing_host_ids = sorted(
                batch_df.loc[batch_df['Host_Sequence'].isna(), 'Host_ID'].drop_duplicates().tolist()
            )

            # Filter out rows with missing sequences
            before_count = len(batch_df)
            batch_df = batch_df.dropna(subset=['Phage_Sequence', 'Host_Sequence'])
            after_count = len(batch_df)

            if before_count > after_count:
                logging.warning(
                    "⚠️  Removed %d pairs with missing sequences from batch "
                    "(%d phages missing, %d hosts missing)",
                    before_count - after_count,
                    len(missing_phage_ids),
                    len(missing_host_ids),
                )
                if missing_phage_ids:
                    logging.warning("   Missing phage IDs (sample): %s", ", ".join(missing_phage_ids[:5]))
                if missing_host_ids:
                    logging.warning("   Missing host IDs (sample): %s", ", ".join(missing_host_ids[:5]))

            if len(batch_df) > 0:
                yield batch_df

        logging.info(f"✅ Completed iteration over {batch_num} batches")
    
    def _get_sequence_safe(
        self,
        seq_id: str,
        seq_type: str,
        host_contig_mode: str = "first",
    ) -> str:
        """Helper method for safe sequence retrieval with error handling.

        Args:
            seq_id: Sequence identifier (Phage_ID or Host_ID).
            seq_type: Type of sequence – ``'phage'`` or ``'host'``.
            host_contig_mode: Contig assembly mode forwarded to
                :meth:`get_host_sequence` when *seq_type* is ``'host'``.
                See :meth:`get_host_sequence` for details.

        Returns:
            Sequence as a string, or an empty string if not found.
        """
        try:
            if seq_type == 'phage':
                return str(self.phage_fasta[seq_id][:].seq)
            elif seq_type == 'host':
                return self.get_host_sequence(seq_id, contig_mode=host_contig_mode)
            else:
                raise ValueError(f"Invalid seq_type: {seq_type}")
        except KeyError:
            logging.warning(f"⚠️  {seq_type.capitalize()} sequence not found for ID: {seq_id}")
            return ""
        except Exception as e:
            logging.warning(f"⚠️  Error retrieving {seq_type} sequence for {seq_id}: {e}")
            return ""
    
    def close(self):
        """Close database connection"""
        self.conn.close()
        logging.info("🔒 Database connection closed")
    
    def _fetch_phage_sequences(self, phage_ids: list) -> pd.DataFrame:
        """
        Fetch phage sequences for given Phage IDs.

        When a private phage mapping is configured, the source database and
        source type are looked up from the DB so that private phages are routed
        to their per-source FASTA rather than all_phages.fasta.
        
        Args:
            phage_ids: List of Phage IDs to retrieve sequences for
            
        Returns:
            DataFrame with Phage_ID and Sequence columns
        """
        if not phage_ids:
            logging.warning("No Phage IDs provided")
            return pd.DataFrame(columns=['Phage_ID', 'Sequence'])
        
        logging.info(f"🔍 Fetching sequences for {len(phage_ids):,} phages")
        
        # Resolve per-phage source information for routing private vs. public lookups.
        phage_source_db: Dict[str, str] = {}
        phage_source_type: Dict[str, str] = {}
        if self._private_phage_mapping and phage_ids:
            try:
                placeholders = ", ".join(["?" for _ in phage_ids])
                source_df = self.conn.execute(
                    f"SELECT Phage_ID, Source_DB, source_type FROM fact_phages "
                    f"WHERE Phage_ID IN ({placeholders})",
                    phage_ids,
                ).fetchdf()
                phage_source_db = dict(zip(source_df["Phage_ID"], source_df["Source_DB"]))
                phage_source_type = dict(zip(source_df["Phage_ID"], source_df["source_type"]))
            except Exception as exc:
                logging.debug("Could not query phage source info: %s", exc)

        # Read sequences, routing private phages to their per-source FASTA.
        sequences = []
        missing_ids = []
        
        for phage_id in phage_ids:
            seq = self._get_phage_sequence(
                phage_id,
                source_db=phage_source_db.get(phage_id),
                source_type=phage_source_type.get(phage_id),
            )
            if seq is not None:
                sequences.append({'Phage_ID': phage_id, 'Sequence': seq})
            else:
                missing_ids.append(phage_id)
                logging.warning(f"⚠️  Phage ID '{phage_id}' not found in FASTA")
        
        if missing_ids:
            logging.warning(f"⚠️  {len(missing_ids):,} phage IDs not found in FASTA file")
        
        df = pd.DataFrame(sequences)
        logging.info(f"✅ Retrieved {len(df):,} sequences")
        
        return df

    def _fetch_protein_sequences(self, protein_ids: list) -> pd.DataFrame:
        """
        Fetch protein sequences for given Protein IDs.
        
        Args:
            protein_ids: List of Protein IDs to retrieve sequences for
            
        Returns:
            DataFrame with Protein_ID and Sequence columns
        """
        if not protein_ids:
            logging.warning("No Protein IDs provided")
            return pd.DataFrame(columns=['Protein_ID', 'Sequence'])
        
        logging.info(f"🔍 Fetching sequences for {len(protein_ids):,} proteins")
        
        sequences = []
        missing_ids = []
        
        for protein_id in protein_ids:
            try:
                seq = self.protein_fasta[protein_id][:].seq
                sequences.append({
                    'Protein_ID': protein_id,
                    'Sequence': str(seq)
                })
            except KeyError:
                # Fallback: the FASTA key is the full header (e.g. "AE002163.1 CDS_1 desc").
                # Try prefix match first, then a token-index lookup.
                seq = self._fuzzy_protein_lookup(protein_id)
                if seq is not None:
                    sequences.append({'Protein_ID': protein_id, 'Sequence': seq})
                else:
                    missing_ids.append(protein_id)
                    logging.warning(f"⚠️  Protein ID '{protein_id}' not found in FASTA")
        
        if missing_ids:
            logging.warning(f"⚠️  {len(missing_ids):,} protein IDs not found in FASTA file")
        
        df = pd.DataFrame(sequences)
        logging.info(f"✅ Retrieved {len(df):,} sequences")
        
        return df

    def _build_protein_token_index(self) -> dict:
        """
        Build (once) a mapping of every whitespace-delimited token in every protein
        FASTA key to that key.  Used by _fuzzy_protein_lookup for O(1) fallback
        lookups instead of O(n) scans.

        Common words (e.g. "hypothetical", "protein") will map to the first key
        that contained them; accession tokens (e.g. "YP_009137915.1") are unique
        and therefore return the correct key.
        """
        if not hasattr(self, '_protein_token_idx'):
            idx: dict = {}
            for key in self.protein_fasta.keys():
                for token in key.split():
                    idx.setdefault(token, key)
            self._protein_token_idx = idx
        return self._protein_token_idx

    def _fuzzy_protein_lookup(self, protein_id: str) -> Optional[str]:
        """
        Try to find a protein sequence whose FASTA header matches protein_id
        approximately.

        Strategy (in order):
        1. Prefix match: the full header starts with protein_id.
        2. Token match: the second whitespace-delimited token of protein_id
           (the protein accession when the first token is the phage accession)
           appears as an exact token in a FASTA header, using a pre-built index
           so the lookup is O(1).

        Returns the sequence string, or None if no match is found.
        """
        fasta = self.protein_fasta

        # 1. Prefix match (handles "AE002163.1 CDS_1" matching key "AE002163.1 CDS_1 desc")
        for key in fasta.keys():
            if key.startswith(protein_id):
                return str(fasta[key][:].seq)

        # 2. Token-index lookup (O(1) after the index is built once)
        pid_parts = protein_id.split()
        # When the first token is the phage accession, use the second token
        # (actual protein accession) as the search target; otherwise use the
        # full protein_id as a single token.
        search_token = pid_parts[1] if len(pid_parts) > 1 else pid_parts[0]
        idx = self._build_protein_token_index()
        matched_key = idx.get(search_token)
        if matched_key is not None:
            return str(fasta[matched_key][:].seq)

        return None

# Example usage and testing
if __name__ == "__main__":
    import sys
    
    # Test script
    if len(sys.argv) < 4:
        print("Usage: sequence_retrieval.py <db_path> <phage_fasta> <protein_fasta>")
        sys.exit(1)
    
    db_path = sys.argv[1]
    phage_fasta = sys.argv[2]
    protein_fasta = sys.argv[3]
    
    print("🧪 Testing SequenceRetriever")
    print("="*80)
    
    try:
        # Initialize
        retriever = SequenceRetriever(db_path, phage_fasta, protein_fasta)
        
        # Get stats
        retriever.get_stats()
        
        # Test 1: Get phages by query
        print("\n🧪 Test 1: Query-based phage retrieval")
        query = "SELECT Phage_ID FROM fact_phages WHERE Length > 50000 LIMIT 5"
        phage_df = retriever.get_phage_sequences(query)
        print(phage_df)
        
        # Test 2: Get proteins by query
        print("\n🧪 Test 2: Query-based protein retrieval")
        query = "SELECT Protein_ID FROM dim_proteins LIMIT 5"
        protein_df = retriever.get_protein_sequences(query)
        print(protein_df)
        
        # Test 3: Export to FASTA
        print("\n🧪 Test 3: Export to FASTA")
        retriever.export_fasta(phage_df, "test_output_phages.fasta", "Phage_ID")
        retriever.export_fasta(protein_df, "test_output_proteins.fasta", "Protein_ID")
        
        print("\n✅ All tests passed!")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    finally:
        retriever.close()
