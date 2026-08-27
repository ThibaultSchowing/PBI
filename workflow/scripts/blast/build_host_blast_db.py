"""Build a BLAST nucleotide database from individual host genome FASTA files."""

import json
import os
import subprocess
import sys

from snakemake.io import touch


def main():
    host_mapping_path = snakemake.input.mapping
    db_dir = snakemake.params.db_dir
    db_name = snakemake.params.db_name
    combined_fasta = snakemake.params.combined_fasta
    output_path = snakemake.output.db
    log_path = snakemake.log[0]

    os.makedirs(db_dir, exist_ok=True)

    # Load host FASTA mapping
    with open(host_mapping_path) as f:
        host_mapping = json.load(f)

    if not host_mapping:
        print("Warning: No host genomes found in mapping. Creating empty BLAST DB.",
              file=sys.stderr)
        # Create an empty FASTA so makeblastdb doesn't fail
        with open(combined_fasta, "w") as f:
            f.write(">empty_placeholder\n")
            f.write("N\n")
    else:
        # Concatenate all host genomes into a single FASTA
        with open(combined_fasta, "w") as out_f:
            for host_id, fasta_path in sorted(host_mapping.items()):
                if not os.path.exists(fasta_path):
                    print(f"Warning: Host FASTA not found: {fasta_path} (host_id={host_id})",
                          file=sys.stderr)
                    continue
                with open(fasta_path) as in_f:
                    for line in in_f:
                        out_f.write(line)

    # Build BLAST database
    db_prefix = os.path.join(db_dir, db_name)
    cmd = [
        "makeblastdb",
        "-in", combined_fasta,
        "-dbtype", "nucl",
        "-out", db_prefix,
        "-parse_seqids",
        "-title", "PBI Host Genomes",
    ]

    with open(log_path, "w") as log_f:
        result = subprocess.run(cmd, stderr=log_f, stdout=log_f)

    if result.returncode != 0:
        print(f"makeblastdb failed with return code {result.returncode}", file=sys.stderr)
        sys.exit(1)

    touch(output_path)


main()
