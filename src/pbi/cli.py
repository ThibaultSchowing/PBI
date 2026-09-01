from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .private_data import validate_private_roots

# Default private data root relative to the current working directory.
_DEFAULT_PRIVATE_ROOT = "private_data"

# Default GFF3 paths (container mode)
_DEFAULT_GFF3_DIR = "/data/processed/gff3"
_DEFAULT_GFF3_INDEX = "/data/processed/gff3/gff3_index.json"


def _cmd_validate_private(args: argparse.Namespace) -> int:
    roots = list(args.path) if args.path else []

    if not roots:
        default = Path(_DEFAULT_PRIVATE_ROOT)
        if default.is_dir():
            roots = [str(default)]
            print(f"No --path given; using default: {default.resolve()}")
        else:
            print(
                f"No --path given and default directory '{_DEFAULT_PRIVATE_ROOT}' not found. "
                "Please run this command from the repository root or pass --path explicitly.",
                file=sys.stderr,
            )
            return 2

    summary = validate_private_roots(roots)

    print("Private dataset validation summary")
    print(f"  Roots: {summary['roots']}")
    print(f"  Sources found: {summary['sources_found']}")
    print(f"  Sources valid: {summary['sources_valid']}")
    print(f"  Sources invalid: {summary['sources_invalid']}")
    print("")

    for source in summary["sources"]:
        status = "VALID" if source.is_valid else "INVALID"
        print(f"[{status}] {source.source_db} ({source.source_dir})")
        if source.warnings:
            for warning in source.warnings:
                print(f"  - warning: {warning}")
        if source.errors:
            for error in source.errors:
                print(f"  - error: {error}")
        else:
            print(
                f"  - rows={source.stats['rows']}, unique_phages={source.stats['unique_phages']}, "
                f"unique_hosts={source.stats['unique_hosts']}"
            )
        print("")

    if summary["sources_found"] == 0:
        print("No private sources found. Add subdirectories under the private data root.")

    return 0 if summary["sources_invalid"] == 0 else 1


# Default BLAST database path (container mode)
_DEFAULT_BLAST_DB_DIR = "/data/processed/blast_db"


def _cmd_blast_search(args: argparse.Namespace) -> int:
    """Search a sequence against BLAST databases."""
    from .blast_search import BlastSearcher

    blast_db_dir = args.blast_db_dir or _DEFAULT_BLAST_DB_DIR

    if not Path(blast_db_dir).exists():
        print(f"BLAST database directory not found: {blast_db_dir}", file=sys.stderr)
        print("Run the pipeline to build BLAST databases first.", file=sys.stderr)
        return 1

    try:
        searcher = BlastSearcher(blast_db_dir)

        # Check databases
        if args.list_databases:
            dbs = searcher.list_databases()
            print("Available BLAST databases:")
            for name, info in dbs.items():
                status = "READY" if info["exists"] else "NOT BUILT"
                print(f"  [{status}] {name} ({info['type']})")
            return 0

        # Validate inputs
        if not args.sequence and not args.input:
            print("Error: provide --sequence or --input", file=sys.stderr)
            return 1

        program = args.program
        db = args.db
        max_hits = args.max_hits
        evalue = args.evalue

        if args.input:
            results = searcher.search_fasta(
                args.input, program=program, db=db,
                max_hits=max_hits, evalue=evalue,
            )
        else:
            results = searcher.search_sequence(
                args.sequence, program=program, db=db,
                max_hits=max_hits, evalue=evalue,
            )

        if results.empty:
            print("No hits found.")
            return 0

        # Output
        if args.output:
            results.to_csv(args.output, index=False)
            print(f"Results written to {args.output} ({len(results)} hits)")
        else:
            print(results.to_string(index=False))

        return 0

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error running BLAST search: {e}", file=sys.stderr)
        return 1


def _cmd_get_gff3(args: argparse.Namespace) -> int:
    """Retrieve GFF3 content for a phage."""
    from .gff3_retrieval import GFF3Retriever

    gff3_dir = args.gff3_dir or _DEFAULT_GFF3_DIR
    index_path = args.index or _DEFAULT_GFF3_INDEX

    if not Path(index_path).exists():
        print(f"GFF3 index not found: {index_path}", file=sys.stderr)
        print("Run the pipeline first to generate the GFF3 index.", file=sys.stderr)
        return 1

    try:
        retriever = GFF3Retriever(gff3_dir, index_path)

        if args.stats:
            stats = retriever.stats()
            print(f"GFF3 Index Statistics:")
            print(f"  Total phages: {stats['total_phages']}")
            print(f"  Sources:")
            for source, count in sorted(stats["sources"].items()):
                print(f"    {source}: {count}")
            return 0

        if args.list_sources:
            stats = retriever.stats()
            for source in sorted(stats["sources"].keys()):
                print(source)
            return 0

        if args.list_phages:
            source_filter = args.list_phages if args.list_phages != "all" else None
            phages = retriever.list_phages(source_db=source_filter)
            for phage_id in sorted(phages):
                print(phage_id)
            return 0

        # Default: retrieve GFF3 content for a phage
        if not args.phage_id:
            print("Error: phage_id is required (or use --stats, --list-sources, --list-phages)", file=sys.stderr)
            return 1

        content = retriever.get_gff3(args.phage_id)
        if content:
            print(content, end="")
            return 0
        else:
            print(f"No GFF3 found for '{args.phage_id}'", file=sys.stderr)
            return 1

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error retrieving GFF3: {e}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pbi", description="PBI command line utilities")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_private = subparsers.add_parser(
        "validate-private",
        help="Validate private datasets without running the full pipeline",
    )
    validate_private.add_argument(
        "--path",
        action="append",
        required=False,
        default=None,
        help=(
            "Private data root directory containing one subdirectory per source. "
            "Can be provided multiple times. "
            f"Defaults to ./{_DEFAULT_PRIVATE_ROOT}/ when run from the repository root."
        ),
    )
    validate_private.set_defaults(func=_cmd_validate_private)

    get_gff3 = subparsers.add_parser(
        "get-gff3",
        help="Retrieve GFF3 annotations for a phage or inspect the GFF3 index",
    )
    get_gff3.add_argument(
        "phage_id",
        nargs="?",
        default=None,
        help="Phage ID to retrieve (e.g., 'Mycobacterium_phage_NuevoMundo')",
    )
    get_gff3.add_argument(
        "--gff3-dir",
        default=None,
        help=f"GFF3 directory path (default: {_DEFAULT_GFF3_DIR})",
    )
    get_gff3.add_argument(
        "--index",
        default=None,
        help=f"GFF3 index path (default: {_DEFAULT_GFF3_INDEX})",
    )
    get_gff3.add_argument(
        "--stats",
        action="store_true",
        help="Show GFF3 index statistics",
    )
    get_gff3.add_argument(
        "--list-sources",
        action="store_true",
        help="List all source databases in the index",
    )
    get_gff3.add_argument(
        "--list-phages",
        nargs="?",
        const="all",
        default=None,
        help="List all phage IDs (optionally filter by source DB)",
    )
    get_gff3.set_defaults(func=_cmd_get_gff3)

    blast_search = subparsers.add_parser(
        "blast-search",
        help="Search a sequence against BLAST databases",
    )
    blast_search.add_argument(
        "sequence",
        nargs="?",
        default=None,
        help="Query sequence (DNA or protein)",
    )
    blast_search.add_argument(
        "--input", "-i",
        default=None,
        help="Path to query FASTA file",
    )
    blast_search.add_argument(
        "--program", "-p",
        default="blastn",
        choices=["blastn", "blastp", "blastx", "tblastn", "tblastx"],
        help="BLAST program (default: blastn)",
    )
    blast_search.add_argument(
        "--db", "-d",
        default=None,
        choices=["phages", "proteins", "hosts"],
        help="Database to search (auto-selected if not given)",
    )
    blast_search.add_argument(
        "--max-hits", "-m",
        type=int,
        default=10,
        help="Maximum number of hits (default: 10)",
    )
    blast_search.add_argument(
        "--evalue", "-e",
        type=float,
        default=1e-5,
        help="E-value threshold (default: 1e-5)",
    )
    blast_search.add_argument(
        "--output", "-o",
        default=None,
        help="Output CSV file (prints to stdout if not given)",
    )
    blast_search.add_argument(
        "--blast-db-dir",
        default=None,
        help="Path to BLAST database directory",
    )
    blast_search.add_argument(
        "--list-databases",
        action="store_true",
        help="List available BLAST databases and their status",
    )
    blast_search.set_defaults(func=_cmd_blast_search)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
