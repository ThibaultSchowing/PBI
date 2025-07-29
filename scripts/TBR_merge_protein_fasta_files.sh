#!/usr/bin/env bash
set -euo pipefail

# Ensure exactly one argument is given: the data directory path
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <data_dir>" >&2
    exit 1
fi

# The top-level data directory (e.g., /path/to/data)
data_dir="$1"
# The directory where protein FASTA files are organized by source
protein_fasta_dir="${data_dir}/protein_fasta"

# Enable glob patterns to return empty arrays instead of literal strings
shopt -s nullglob

# Loop over all subdirectories inside protein_fasta_dir (i.e., the sources)
for source_path in "$protein_fasta_dir"/*; do
    # Skip if it's not a directory (e.g., ignore files)
    [[ -d "$source_path" ]] || continue

    # Extract the source name (e.g., "source1" from ".../protein_fasta/source1")
    source=$(basename "$source_path")
    # Expected nested structure: protein_fasta/source/source/phage/
    nested_dir="$source_path/$source"

    # Initialize the list of FASTA files
    fasta_files=()

    # If the nested directory exists, collect all .fasta or .fa files within it
    if [[ -d "$nested_dir" ]]; then
        # Read all matching files recursively from nested_dir
        while IFS= read -r -d '' file; do
            fasta_files+=("$file")
        done < <(find "$nested_dir" -type f \( -iname "*.fasta" -o -iname "*.fa" \) -print0)
    else
        # If no nested directory, try flat layout (directly in source_path)
        fasta_files=("$source_path"/*.fa "$source_path"/*.fasta)
    fi

    # Filter out non-existing entries (glob patterns might yield "*.fasta")
    real_files=()
    for f in "${fasta_files[@]}"; do
        [[ -f "$f" ]] && real_files+=("$f")
    done

    file_count="${#real_files[@]}"

    # Handle cases based on the number of FASTA files found
    if [[ "$file_count" -eq 0 ]]; then
        echo "⚠️  No FASTA files found for $source" >&2
        continue

    elif [[ "$file_count" -eq 1 ]]; then
        echo "✅ Single FASTA file for $source — no merge needed: ${real_files[0]}"
        continue

    else
        # Define the output merged file path: protein_fasta/source/source_merged.fasta
        merged="${source_path}/${source}_merged.fasta"
        > "$merged"  # Truncate or create the merged file

        echo "🔄 Merging $file_count FASTA files for $source → $merged"

        # Loop through each input file and copy with modified headers
        for f in "${real_files[@]}"; do
            # The phage name is assumed to be the parent directory of the file
            phage=$(basename "$(dirname "$f")")

            # Use awk to modify headers by prepending phage name
            awk -v phage="$phage" '
                /^>/ { sub(/^>/, ">" phage " "); print; next }  # Modify headers
                { print }                                       # Keep sequences unchanged
            ' "$f" >> "$merged"
        done

        echo "✅ Done merging for $source"
    fi
done

