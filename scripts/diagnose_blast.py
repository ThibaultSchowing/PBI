#!/usr/bin/env python3
"""
Diagnostic script for BLAST database issues.

Usage:
    python scripts/diagnose_blast.py
    python scripts/diagnose_blast.py --db phages
    python scripts/diagnose_blast.py --test-search
    python scripts/diagnose_blast.py --db phages --test-search --timeout 600

In Docker:
    docker exec pbi-analysis python scripts/diagnose_blast.py
    docker exec pbi-analysis python scripts/diagnose_blast.py --db phages --test-search
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pbi import BlastSearcher


def format_size(size_bytes: int) -> str:
    """Format size in human-readable format."""
    if size_bytes >= 1024 ** 3:
        return f"{size_bytes / (1024 ** 3):.2f} GB"
    elif size_bytes >= 1024 ** 2:
        return f"{size_bytes / (1024 ** 2):.1f} MB"
    elif size_bytes >= 1024:
        return f"{size_bytes / 1024:.1f} KB"
    else:
        return f"{size_bytes} bytes"


def diagnose_database(searcher: BlastSearcher, db_name: str, test_search: bool = False,
                      timeout: int = 600, num_threads: int = 1):
    """Diagnose a specific BLAST database."""
    print(f"\n{'='*60}")
    print(f"Database: {db_name}")
    print(f"{'='*60}")

    # Get database info from list_databases
    dbs = searcher.list_databases()
    if db_name not in dbs:
        print(f"ERROR: Database '{db_name}' not found in list_databases()")
        return

    db_info = dbs[db_name]
    print(f"Type: {db_info['type']}")
    print(f"Exists: {db_info['exists']}")
    print(f"Path: {db_info['path']}")
    print(f"Total Size: {format_size(int(db_info.get('total_size_mb', 0) * 1024 * 1024))}")

    # Test search if requested
    if test_search:
        print(f"\nTesting BLAST search...")
        test_query = "GTTCTTGTCGAAAAACGTCAACATTTTATAAAAAAGGGTTGCA"
        print(f"Query: {test_query}")
        print(f"Query length: {len(test_query)} bp")
        print(f"Timeout: {timeout} seconds")
        print(f"Threads: {num_threads}")

        start_time = time.time()
        try:
            results = searcher.search_sequence(
                test_query,
                program="blastn",
                db=db_name,
                max_hits=5,
                evalue=1e-5,
                timeout=timeout,
                num_threads=num_threads,
            )
            elapsed = time.time() - start_time
            print(f"SUCCESS: Search completed in {elapsed:.2f} seconds")
            print(f"Found {len(results)} hits")
            if not results.empty:
                print(f"Top hit: {results.iloc[0]['sseqid']} ({results.iloc[0]['pident']}%)")
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"FAILED: Search failed after {elapsed:.2f} seconds")
            print(f"Error: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose BLAST database issues",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Diagnose all databases
  python scripts/diagnose_blast.py

  # Diagnose specific database
  python scripts/diagnose_blast.py --db phages

  # Test search with extended timeout
  python scripts/diagnose_blast.py --db phages --test-search

  # Test search with custom timeout and threads
  python scripts/diagnose_blast.py --db phages --test-search --timeout 1200 --num-threads 4

In Docker:
  docker exec pbi-analysis python scripts/diagnose_blast.py --db phages --test-search
        """
    )
    parser.add_argument(
        "--db",
        default=None,
        choices=["phages", "proteins", "hosts", "private", "combined"],
        help="Specific database to diagnose (default: all)"
    )
    parser.add_argument(
        "--test-search",
        action="store_true",
        help="Test BLAST search with a sample query"
    )
    parser.add_argument(
        "--timeout", "-t",
        type=int,
        default=600,
        help="Timeout for test search in seconds (default: 600 = 10 minutes)"
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=1,
        help="Number of threads for test search (default: 1)"
    )
    parser.add_argument(
        "--blast-db-dir",
        default=None,
        help="Path to BLAST database directory"
    )
    args = parser.parse_args()

    # Initialize searcher
    if args.blast_db_dir:
        searcher = BlastSearcher(args.blast_db_dir)
    else:
        searcher = BlastSearcher()

    print("BLAST Diagnostic Tool")
    print("=" * 60)
    print(f"BLAST binary directory: {searcher.blast_bin_dir}")
    print(f"BLAST database directory: {searcher.blast_db_dir}")

    # List all databases
    print(f"\nAll Databases:")
    print("-" * 60)

    dbs = searcher.list_databases()
    for name, info in dbs.items():
        status = "READY" if info["exists"] else "NOT BUILT"
        size_mb = info.get("total_size_mb", 0)
        size_gb = info.get("total_size_gb", 0)
        if size_gb >= 1:
            size_str = f"{size_gb:.2f} GB"
        else:
            size_str = f"{size_mb:.1f} MB"
        print(f"  [{status}] {name} ({info['type']}) - {size_str}")

    # Diagnose specific database(s)
    if args.db:
        diagnose_database(searcher, args.db, args.test_search, args.timeout, args.num_threads)
    else:
        # Diagnose all databases (without test search unless requested)
        for name in dbs.keys():
            if dbs[name]["exists"]:
                diagnose_database(searcher, name, args.test_search, args.timeout, args.num_threads)

    print(f"\n{'='*60}")
    print("Diagnostic complete")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
