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
            type (nucl/prot), exists (bool), path (str).
        """
        databases = {}
        for db_name, db_type in [
            ("phages", "nucl"),
            ("proteins", "prot"),
            ("hosts", "nucl"),
        ]:
            db_path = self.blast_db_dir / db_name
            done_marker = db_path / f"makeblastdb_{db_name}.done"
            db_prefix = db_path / f"all_{db_name}"
            exists = done_marker.exists() or (
                db_prefix.with_suffix(".nsq").exists()  # nucleotide
                if db_type == "nucl"
                else db_prefix.with_suffix(".psq").exists()  # protein
            )
            databases[db_name] = {
                "type": db_type,
                "exists": exists,
                "path": str(db_path),
            }
        return databases

    def get_db_prefix(self, db_name: str) -> str:
        """
        Get the BLAST database prefix for a named database.

        Args:
            db_name: One of 'phages', 'proteins', 'hosts'.

        Returns:
            Full path prefix for the BLAST database.
        """
        valid_dbs = {"phages", "proteins", "hosts"}
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

    def search_sequence(
        self,
        sequence: str,
        program: str = "blastn",
        db: str | None = None,
        max_hits: int = 10,
        evalue: float = 1e-5,
        outfmt: int = 6,
        extra_args: list[str] | None = None,
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

        return self._run_blast(
            query=str(fasta_path),
            program=program,
            db=db,
            max_hits=max_hits,
            evalue=evalue,
            outfmt=outfmt,
            extra_args=extra_args,
        )

    def _run_blast(
        self,
        query: str,
        program: str,
        db: str,
        max_hits: int,
        evalue: float,
        outfmt: int,
        extra_args: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Execute a BLAST search and return results as a DataFrame.
        """
        db_prefix = self.get_db_prefix(db)

        columns = BLAST_OUTFMT7_COLUMNS if outfmt == 7 else BLAST_OUTFMT6_COLUMNS
        fmt_str = f"{outfmt} {' '.join(columns)}"

        cmd = [
            self._blast_bin(program),
            "-query", query,
            "-db", db_prefix,
            "-outfmt", fmt_str,
            "-max_target_seqs", str(max_hits),
            "-evalue", str(evalue),
            "-num_threads", "1",
        ]

        if extra_args:
            cmd.extend(extra_args)

        logger.debug("Running BLAST: %s", " ".join(cmd))

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
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

        Returns:
            DataFrame with BLAST results and annotations.
        """
        df = self.search_sequence(
            sequence, program=program, db=db,
            max_hits=max_hits, evalue=evalue, outfmt=7,
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
