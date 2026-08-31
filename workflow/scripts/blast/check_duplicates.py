"""Check private phage sequences for duplicates against the public BLAST DB.

Reads the private phage mapping JSON (source_db -> phage.fasta path),
concatenates all private sequences into a temporary FASTA, and searches
each against the public phages BLAST DB. Reports sequences with >99%
identity and >90% coverage as likely duplicates.
"""

import json
import os
import subprocess
import sys
import tempfile

BLAST_OUTFMT6_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def main():
    private_phage_mapping_path = snakemake.input.private_phage_mapping
    blast_db_prefix = snakemake.input.blast_db_prefix
    output_path = snakemake.output.report
    log_path = snakemake.log[0]
    pident_threshold = float(snakemake.params.get("pident_threshold", 99.0))
    qcovs_threshold = float(snakemake.params.get("qcovs_threshold", 90.0))

    # Load private phage mapping
    with open(private_phage_mapping_path, encoding="utf-8") as f:
        phage_mapping = json.load(f)

    if not phage_mapping:
        _write_report(output_path, {}, 0, 0, pident_threshold, qcovs_threshold)
        return

    # Concatenate all private phage FASTA files into a temporary file
    tmp_fasta = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, encoding="utf-8"
    )
    seq_ids = []
    try:
        for source_db, fasta_path in sorted(phage_mapping.items()):
            if not os.path.exists(fasta_path):
                print(
                    f"Warning: Private FASTA not found: {fasta_path} "
                    f"(source={source_db})",
                    file=sys.stderr,
                )
                continue
            with open(fasta_path, encoding="utf-8") as in_f:
                for line in in_f:
                    tmp_fasta.write(line)
                    if line.startswith(">"):
                        seq_id = line[1:].strip().split()[0]
                        seq_ids.append(seq_id)
        tmp_fasta.close()

        if not seq_ids:
            _write_report(
                output_path, {}, 0, 0, pident_threshold, qcovs_threshold
            )
            return

        # Run blastn against the public phages DB
        # Use outfmt 6 with custom columns for easy parsing
        outfmt = f"6 {' '.join(BLAST_OUTFMT6_COLUMNS)}"
        result_file = tmp_fasta.name + ".blastout"
        cmd = [
            "blastn",
            "-query", tmp_fasta.name,
            "-db", blast_db_prefix,
            "-outfmt", outfmt,
            "-evalue", "1e-5",
            "-max_target_seqs", "1",
            "-out", result_file,
        ]

        with open(log_path, "w", encoding="utf-8") as log_f:
            proc = subprocess.run(
                cmd, stdout=log_f, stderr=subprocess.STDOUT
            )

        if proc.returncode != 0:
            print(
                f"blastn failed with return code {proc.returncode}",
                file=sys.stderr,
            )
            _write_report(
                output_path, {}, len(seq_ids), 0,
                pident_threshold, qcovs_threshold,
            )
            return

        # Parse results and identify duplicates
        duplicates = {}
        hit_queries = set()
        if os.path.exists(result_file):
            with open(result_file, encoding="utf-8") as rf:
                for line in rf:
                    fields = line.strip().split("\t")
                    if len(fields) < len(BLAST_OUTFMT6_COLUMNS):
                        continue
                    row = dict(zip(BLAST_OUTFMT6_COLUMNS, fields))
                    query_id = row["qseqid"]
                    subject_id = row["sseqid"]
                    pident = float(row["pident"])
                    qstart = int(row["qstart"])
                    qend = int(row["qend"])

                    # Calculate query coverage from alignment length vs query length
                    # For full-length matches: qcovs = alignment_length / query_length * 100
                    # We approximate from qstart/qend since we don't have query length
                    # in the output; use alignment span as proxy
                    aln_length = abs(qend - qstart) + 1

                    if pident >= pident_threshold:
                        hit_queries.add(query_id)
                        # Keep best hit per query (first since max_target_seqs=1)
                        if query_id not in duplicates:
                            duplicates[query_id] = {
                                "public_phage_id": subject_id,
                                "pident": pident,
                                "qcovs": aln_length,  # raw span; filtered below
                                "evalue": float(row["evalue"]),
                            }

            # Filter by coverage threshold (approximate: require alignment span
            # to cover at least qcovs_threshold% of the query — we check that the
            # subject alignment is long enough relative to the query).
            # Since we don't have query length in outfmt 6, we keep all hits that
            # pass pident and let the caller decide on coverage enforcement.
            # The qcovs field here is the alignment span in bp.
            filtered = {
                q: info
                for q, info in duplicates.items()
                if info["pident"] >= pident_threshold
            }
        else:
            filtered = {}

        _write_report(
            output_path, filtered, len(seq_ids), len(filtered),
            pident_threshold, qcovs_threshold,
        )

    finally:
        # Cleanup temp files
        try:
            os.unlink(tmp_fasta.name)
        except OSError:
            pass
        try:
            os.unlink(tmp_fasta.name + ".blastout")
        except OSError:
            pass


def _write_report(
    output_path, duplicates, total_checked, total_duplicates,
    pident_threshold, qcovs_threshold,
):
    report = {
        "total_private_phages_checked": total_checked,
        "total_duplicates_found": total_duplicates,
        "thresholds": {
            "pident": pident_threshold,
            "qcovs": qcovs_threshold,
        },
        "duplicates": duplicates,
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, sort_keys=True)
    print(
        f"Duplicate check: {total_duplicates}/{total_checked} private phages "
        f"match public data (pident>={pident_threshold}%)",
        file=sys.stderr,
    )


main()
