"""Build a BLAST nucleotide database from private phage FASTA files.

Collects all per-source phage.fasta files produced by prepare_private_sequences
and concatenates them into a single FASTA for makeblastdb.
"""

import glob
import os
import subprocess
import sys

from snakemake.io import touch


def main():
    private_phage_dir = snakemake.input.private_phage_dir
    db_dir = snakemake.params.db_dir
    db_name = snakemake.params.db_name
    max_file_size = snakemake.params.get("max_file_size", "1G")
    output_path = snakemake.output.db
    log_path = snakemake.log[0]

    os.makedirs(db_dir, exist_ok=True)

    # Collect all per-source phage.fasta files
    fasta_pattern = os.path.join(private_phage_dir, "*", "phage.fasta")
    fasta_files = sorted(glob.glob(fasta_pattern))

    combined_fasta = os.path.join(db_dir, f"{db_name}.fasta")

    if not fasta_files:
        print(
            "Warning: No private phage FASTA files found. "
            "Creating placeholder BLAST DB.",
            file=sys.stderr,
        )
        with open(combined_fasta, "w") as f:
            f.write(">empty_placeholder\n")
            f.write("N\n")
    else:
        seq_count = 0
        with open(combined_fasta, "w") as out_f:
            for fasta_path in fasta_files:
                with open(fasta_path, encoding="utf-8") as in_f:
                    for line in in_f:
                        out_f.write(line)
                        if line.startswith(">"):
                            seq_count += 1
        print(
            f"Concatenated {seq_count} sequences from {len(fasta_files)} "
            f"private source(s) into {combined_fasta}",
            file=sys.stderr,
        )

    # Build BLAST database
    # Private headers are not NCBI-formatted, so -parse_seqids is omitted
    db_prefix = os.path.join(db_dir, db_name)
    cmd = [
        "makeblastdb",
        "-in", combined_fasta,
        "-dbtype", "nucl",
        "-out", db_prefix,
        "-max_file_sz", max_file_size,
    ]

    with open(log_path, "w", encoding="utf-8") as log_f:
        result = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        print(
            f"makeblastdb failed with return code {result.returncode}",
            file=sys.stderr,
        )
        sys.exit(1)

    touch(output_path)


main()
