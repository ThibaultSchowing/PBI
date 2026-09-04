"""
BLAST search wrapper for PBI.

Provides a Python interface to NCBI BLAST+ for searching unknown sequences
against pre-built BLAST databases (phage genomes, proteins, host genomes).
"""

from __future__ import annotations

import hashlib
import io
import logging
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# BLAST tabular output format 6 columns
BLAST_OUTFMT6_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

# Extended format 7 columns (with query/subject length)
BLAST_OUTFMT7_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovs", "scovhsp",
]

# Default BLAST program to database mapping
DEFAULT_DB_MAP = {
    "blastn": "phages",
    "blastp": "proteins",
    "blastx": "proteins",
    "tblastn": "phages",
    "tblastx": "phages",
}

# Valid BLAST programs
VALID_PROGRAMS = {"blastn", "blastp", "blastx", "tblastn", "tblastx"}


def _find_blast_bin() -> str:
    """Find the BLAST+ executable directory."""
    blastn = shutil.which("blastn")
    if blastn:
        return str(Path(blastn).parent)
    # Fallback: common conda/system locations
    for candidate in ["/usr/bin", "/usr/local/bin", "/opt/conda/bin"]:
        if os.path.isfile(os.path.join(candidate, "blastn")):
            return candidate
    raise FileNotFoundError(
        "BLAST+ executables not found. Install via conda: conda install -c bioconda blast"
    )


class BlastSearcher:
    """
    Wrapper for NCBI BLAST+ searches against PBI databases.

    Usage::

        searcher = BlastSearcher("/data/processed/blast_db")
        results = searcher.search_sequence("ATGCGTTTACG...", program="blastn", db="phages")
    """

    def __init__(self, blast_db_dir: str | Path | None = None):
        """
        Initialize BlastSearcher.

        Args:
            blast_db_dir: Path to the directory containing BLAST databases.
                          Defaults to /data/processed/blast_db (Docker) or
                          <project>/data/processed/blast_db (local).
        """
        if blast_db_dir is None:
            data_path = os.environ.get("DATA_PATH")
            if data_path:
                blast_db_dir = Path(data_path) / "blast_db"
            else:
                project_root = Path(__file__).parent.parent.parent
                blast_db_dir = project_root / "data" / "processed" / "blast_db"

        self.blast_db_dir = Path(blast_db_dir)
        self._blast_bin_dir = None

    @property
    def blast_bin_dir(self) -> str:
        if self._blast_bin_dir is None:
            self._blast_bin_dir = _find_blast_bin()
        return self._blast_bin_dir

    def _blast_bin(self, program: str) -> str:
        """Get the full path to a BLAST+ executable."""
        return os.path.join(self.blast_bin_dir, program)

    def list_databases(self) -> dict[str, dict]:
        """
        List available BLAST databases.

        Returns:
            Dict mapping database name to info dict with keys:
            type (nucl/prot), exists (bool), path (str),
            total_size_mb (float), total_size_gb (float).
        """
        databases = {}
        for db_name, db_type in [
            ("phages", "nucl"),
            ("proteins", "prot"),
            ("hosts", "nucl"),
            ("private", "nucl"),
            ("combined", "nucl"),
        ]:
            db_path = self.blast_db_dir / db_name
            done_marker = db_path / f"makeblastdb_{db_name}.done"
            db_prefix = db_path / f"all_{db_name}"

            # Check for split volumes (e.g., all_phages.00.nsq, all_phages.01.nsq)
            has_volumes = any(
                db_prefix.with_suffix(f".{i:02d}.nsq").exists()
                if db_type == "nucl"
                else db_prefix.with_suffix(f".{i:02d}.psq").exists()
                for i in range(100)
            )

            exists = done_marker.exists() or has_volumes or (
                db_prefix.with_suffix(".nsq").exists()  # nucleotide
                if db_type == "nucl"
                else db_prefix.with_suffix(".psq").exists()  # protein
            )

            total_size = self.get_database_size(db_name)

            databases[db_name] = {
                "type": db_type,
                "exists": exists,
                "path": str(db_path),
                "total_size_mb": total_size / (1024 * 1024),
                "total_size_gb": total_size / (1024 ** 3),
            }
        return databases

    def get_db_prefix(self, db_name: str) -> str:
        """
        Get the BLAST database prefix for a named database.

        Args:
            db_name: One of 'phages', 'proteins', 'hosts', 'private', 'combined'.

        Returns:
            Full path prefix for the BLAST database.
        """
        valid_dbs = {"phages", "proteins", "hosts", "private", "combined"}
        if db_name not in valid_dbs:
            raise ValueError(f"Invalid database '{db_name}'. Must be one of: {valid_dbs}")

        db_path = self.blast_db_dir / db_name / f"all_{db_name}"
        if not db_path.exists():
            # Check if done marker exists (DB was built)
            done_marker = self.blast_db_dir / db_name / f"makeblastdb_{db_name}.done"
            if not done_marker.exists():
                raise FileNotFoundError(
                    f"BLAST database '{db_name}' not found at {db_path}. "
                    "Run the pipeline to build BLAST databases first."
                )
        return str(db_path)

    def get_database_size(self, db_name: str) -> int:
        """
        Get total size of a BLAST database in bytes.

        Handles both single-file and split (volumed) databases.

        Args:
            db_name: One of 'phages', 'proteins', 'hosts', 'private', 'combined'.

        Returns:
            Total size in bytes.
        """
        db_path = self.blast_db_dir / db_name
        db_prefix = db_path / f"all_{db_name}"

        total_size = 0

        # Check for split volumes (.00, .01, .02, etc.)
        for vol_num in range(100):  # Support up to 100 volumes
            vol_suffix = f".{vol_num:02d}"
            for ext in [".nsq", ".psq", ".nhr", ".nin", ".phr", ".pin"]:
                vol_file = db_prefix.with_suffix(vol_suffix + ext)
                if vol_file.exists():
                    total_size += vol_file.stat().st_size

        # If no split volumes found, check single file
        if total_size == 0:
            for ext in [".nsq", ".psq"]:
                single_file = db_prefix.with_suffix(ext)
                if single_file.exists():
                    total_size += single_file.stat().st_size

        return total_size

    def search_sequence(
        self,
        sequence: str,
        program: str = "blastn",
        db: str | None = None,
        max_hits: int = 10,
        evalue: float = 1e-5,
        outfmt: int = 6,
        extra_args: list[str] | None = None,
        timeout: int = 1800,
        num_threads: int = 1,
    ) -> pd.DataFrame:
        """
        Search a single sequence against a BLAST database.

        Args:
            sequence: Query sequence (DNA or protein).
            program: BLAST program (blastn, blastp, blastx, tblastn, tblastx).
            db: Database name (phages, proteins, hosts). Auto-selected if None.
            max_hits: Maximum number of hits to return.
            evalue: E-value threshold.
            outfmt: Output format (6 or 7). Format 7 includes query/subject lengths.
            extra_args: Additional command-line arguments.
            timeout: Timeout in seconds (default 1800 = 30 minutes).
            num_threads: Number of threads for BLAST (default 1).

        Returns:
            DataFrame with BLAST results.
        """
        if program not in VALID_PROGRAMS:
            raise ValueError(f"Invalid program '{program}'. Must be one of: {VALID_PROGRAMS}")

        if db is None:
            db = DEFAULT_DB_MAP.get(program, "phages")

        # Write sequence to temp FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            tmp.write(">query\n")
            tmp.write(sequence.strip() + "\n")
            tmp_path = tmp.name

        try:
            return self._run_blast(
                query=tmp_path,
                program=program,
                db=db,
                max_hits=max_hits,
                evalue=evalue,
                outfmt=outfmt,
                extra_args=extra_args,
                timeout=timeout,
                num_threads=num_threads,
            )
        finally:
            os.unlink(tmp_path)

    def search_fasta(
        self,
        fasta_path: str | Path,
        program: str = "blastn",
        db: str | None = None,
        max_hits: int = 10,
        evalue: float = 1e-5,
        outfmt: int = 6,
        extra_args: list[str] | None = None,
        timeout: int = 1800,
        num_threads: int = 1,
    ) -> pd.DataFrame:
        """
        Search a FASTA file against a BLAST database.

        Args:
            fasta_path: Path to query FASTA file.
            program: BLAST program.
            db: Database name. Auto-selected if None.
            max_hits: Maximum hits per query sequence.
            evalue: E-value threshold.
            outfmt: Output format (6 or 7).
            extra_args: Additional command-line arguments.
            timeout: Timeout in seconds (default 1800 = 30 minutes).
            num_threads: Number of threads for BLAST (default 1).

        Returns:
            DataFrame with BLAST results.
        """
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"Query FASTA not found: {fasta_path}")

        if program not in VALID_PROGRAMS:
            raise ValueError(f"Invalid program '{program}'. Must be one of: {VALID_PROGRAMS}")

        if db is None:
            db = DEFAULT_DB_MAP.get(program, "phages")

        # Normalize line endings — BLAST on Linux cannot parse \r\n
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            content = fasta_path.read_text()
            tmp.write(content.replace("\r\n", "\n").replace("\r", "\n"))
            tmp_path = tmp.name

        try:
            return self._run_blast(
                query=tmp_path,
                program=program,
                db=db,
                max_hits=max_hits,
                evalue=evalue,
                outfmt=outfmt,
                extra_args=extra_args,
                timeout=timeout,
                num_threads=num_threads,
            )
        finally:
            os.unlink(tmp_path)

    def _run_blast(
        self,
        query: str,
        program: str,
        db: str,
        max_hits: int,
        evalue: float,
        outfmt: int,
        extra_args: list[str] | None = None,
        timeout: int = 1800,
        num_threads: int = 1,
    ) -> pd.DataFrame:
        """
        Execute a BLAST search and return results as a DataFrame.

        Args:
            query: Path to query FASTA file.
            program: BLAST program (blastn, blastp, blastx, tblastn, tblastx).
            db: Database name (phages, proteins, hosts, private, combined).
            max_hits: Maximum number of hits to return.
            evalue: E-value threshold.
            outfmt: Output format (6 or 7).
            extra_args: Additional command-line arguments.
            timeout: Timeout in seconds (default 1800 = 30 minutes).
            num_threads: Number of threads for BLAST (default 1).

        Returns:
            DataFrame with BLAST results.
        """
        db_prefix = self.get_db_prefix(db)
        db_size = self.get_database_size(db)
        db_size_gb = db_size / (1024 ** 3)

        logger.info(
            "Starting BLAST search: program=%s, db=%s (%.2f GB), "
            "max_hits=%d, evalue=%g, timeout=%ds, num_threads=%d",
            program, db, db_size_gb, max_hits, evalue, timeout, num_threads,
        )

        columns = BLAST_OUTFMT7_COLUMNS if outfmt == 7 else BLAST_OUTFMT6_COLUMNS
        fmt_str = f"{outfmt} {' '.join(columns)}"

        cmd = [
            self._blast_bin(program),
            "-query", query,
            "-db", db_prefix,
            "-outfmt", fmt_str,
            "-max_target_seqs", str(max_hits),
            "-evalue", str(evalue),
            "-num_threads", str(num_threads),
        ]

        if extra_args:
            cmd.extend(extra_args)

        logger.debug("Running BLAST: %s", " ".join(cmd))

        start_time = time.time()

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as e:
            elapsed = time.time() - start_time
            logger.error(
                "BLAST timed out after %.1f seconds (timeout=%ds). "
                "Database: %s (%.2f GB), Program: %s",
                elapsed, timeout, db, db_size_gb, program,
            )
            raise RuntimeError(
                f"BLAST search timed out after {elapsed:.1f} seconds (timeout={timeout}s).\n"
                f"Database: {db} ({db_size_gb:.2f} GB)\n"
                f"Program: {program}\n"
                f"Query: {query}\n\n"
                f"This database is very large ({db_size_gb:.2f} GB). "
                f"Consider:\n"
                f"1. Increasing timeout: searcher.search_sequence(..., timeout=3600)\n"
                f"2. Using a smaller database\n"
                f"3. Using more threads: searcher.search_sequence(..., num_threads=4)\n"
            ) from e

        elapsed = time.time() - start_time
        logger.info(
            "BLAST completed in %.2f seconds (database: %s, %.2f GB)",
            elapsed, db, db_size_gb,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"BLAST failed (exit code {result.returncode}):\n"
                f"stderr: {result.stderr}\n"
                f"stdout: {result.stdout}"
            )

        # Parse tabular output
        if not result.stdout.strip():
            return pd.DataFrame(columns=columns)

        df = pd.read_csv(
            io.StringIO(result.stdout.strip()),
            sep="\t",
            header=None,
            names=columns,
        )

        return df

    def search_and_annotate(
        self,
        sequence: str,
        program: str = "blastn",
        db: str | None = None,
        max_hits: int = 10,
        evalue: float = 1e-5,
        timeout: int = 1800,
        num_threads: int = 1,
    ) -> pd.DataFrame:
        """
        Search a sequence and add human-readable annotations from the database.

        This enriches BLAST results with metadata from the PBI database
        by extracting the subject sequence identifier.

        Args:
            sequence: Query sequence.
            program: BLAST program.
            db: Database name.
            max_hits: Maximum hits.
            evalue: E-value threshold.
            timeout: Timeout in seconds (default 1800 = 30 minutes).
            num_threads: Number of threads for BLAST (default 1).

        Returns:
            DataFrame with BLAST results and annotations.
        """
        df = self.search_sequence(
            sequence, program=program, db=db,
            max_hits=max_hits, evalue=evalue, outfmt=7,
            timeout=timeout, num_threads=num_threads,
        )

        if df.empty:
            return df

        # Extract subject ID prefix as a source indicator
        df["source_db"] = df["sseqid"].apply(_extract_source_db)

        return df

    @staticmethod
    def query_hash(sequence: str, program: str, db: str, evalue: float) -> str:
        """Generate a hash key for caching BLAST queries."""
        key = f"{sequence.strip().upper()}|{program}|{db}|{evalue}"
        return hashlib.sha256(key.encode()).hexdigest()[:16]


def _extract_source_db(subject_id: str) -> str:
    """Extract likely source database from a sequence identifier."""
    subject_id = str(subject_id)
    known_sources = [
        "RefSeq", "GenBank", "GB", "NCBI", "PhagesDB", "GPD", "GVD",
        "MGV", "TemPhD", "CHVD", "IGVD", "IMGVR", "GOV2", "STV",
        "GSV", "UHGV", "HPGC", "URPC", "ELGV", "OVD", "VMGC",
        "OPD", "BAPS", "SVD", "TYMEFLIES", "MetaVR",
    ]
    for source in known_sources:
        if source.lower() in subject_id.lower():
            return source
    return "Unknown"
