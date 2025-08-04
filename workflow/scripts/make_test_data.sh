#!/usr/bin/env bash
set -euo pipefail

# Base test directory
TEST_BASE="./test_data"
PROTEIN_FASTA="${TEST_BASE}/protein_fasta"

# Remove old test data
rm -rf "$TEST_BASE"
mkdir -p "$PROTEIN_FASTA"

# Function to create a source with phage subdirs and fasta files
create_source() {
    local source="$1"
    local phage_names=("${@:2}")  # All remaining args are phage names
    local source_dir="${PROTEIN_FASTA}/${source}/${source}"

    mkdir -p "$source_dir"

    for phage in "${phage_names[@]}"; do
        local phage_dir="${source_dir}/${phage}"
        mkdir -p "$phage_dir"
        for i in {1..2}; do
            echo -e ">${phage}_seq${i}\nATGCATGCATGC" > "${phage_dir}/${phage}_${i}.fasta"
        done
    done
}

# Create 3 sources:
# - source1 with 2 phages (4 fasta files)
# - source2 with 1 fasta file only (flat layout)
# - source3 with 0 fasta files

create_source "source1" "phageA" "phageB"

# source2: single flat file
mkdir -p "${PROTEIN_FASTA}/source2"
echo -e ">single_seq\nGGTTAACC" > "${PROTEIN_FASTA}/source2/only.fasta"

# source3: empty
mkdir -p "${PROTEIN_FASTA}/source3/source3"

echo "✅ Test directory created at: $TEST_BASE"
tree "$TEST_BASE"

