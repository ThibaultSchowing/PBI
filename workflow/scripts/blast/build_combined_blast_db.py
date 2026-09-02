"""Build a combined BLAST nucleotide database from public + private phage FASTAs.

Merges all_phages.fasta (public) with all per-source private phage.fasta files
into a single FASTA, then runs makeblastdb. This DB is intended for user queries
via the API/CLI so they can search against everything at once.
"""

import glob
import os
import subprocess
import sys

from snakemake.io import touch


def main():
    public_fasta = snakemake.input.public_fasta
    private_phage_dir = snakemake.params.private_phage_dir
    db_dir = snakemake.params.db_dir
    db_name = snakemake.params.db_name
    max_file_size = snakemake.params.get("max_file_size", "1G")
    output_path = snakemake.output.db
    log_path = snakemake.log[0]

    os.makedirs(db_dir, exist_ok=True)

    combined_fasta = os.path.join(db_dir, f"{db_name}.fasta")
    public_count = 0
    private_count = 0

    # Start with public sequences
    with open(combined_fasta, "w", encoding="utf-8") as out_f:
        if os.path.exists(public_fasta):
            with open(public_fasta, encoding="utf-8") as in_f:
                for line in in_f:
                    out_f.write(line)
                    if line.startswith(">"):
                        public_count += 1
        else:
            print(
                f"Warning: Public FASTA not found at {public_fasta}",
                file=sys.stderr,
            )

        # Append private sequences
        fasta_pattern = os.path.join(private_phage_dir, "*", "phage.fasta")
        fasta_files = sorted(glob.glob(fasta_pattern))
        for fasta_path in fasta_files:
            with open(fasta_path, encoding="utf-8") as in_f:
                for line in in_f:
                    out_f.write(line)
                    if line.startswith(">"):
                        private_count += 1

    print(
        f"Combined DB: {public_count} public + {private_count} private = "
        f"{public_count + private_count} total sequences",
        file=sys.stderr,
    )

    # Build BLAST database
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
