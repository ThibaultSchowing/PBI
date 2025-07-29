#!.pixi/envs/default/bin/python

import sys
from pathlib import Path

def merge_fasta_files(source_dir: Path, output_file: Path):
    source_dir = Path(source_dir)
    with output_file.open('w') as out_f:
        for fasta_file in sorted(source_dir.rglob("*.fa")) + sorted(source_dir.rglob("*.fasta")):
            if not fasta_file.is_file():
                continue
            #phage = fasta_file.parent.name # copied from protein script, but not needed for phage fasta as no phage name is added to the header
            
            with fasta_file.open('r') as in_f:
                for line in in_f:
                    if line.startswith('>'):
                        out_f.write(f">{line[1:]}") # shorter than protein as no phage name is needed (already in the sequence header)
                    else:
                        out_f.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: merge_fastas.py <source_dir> <output_file>", file=sys.stderr)
        sys.exit(1)

    src_dir = Path(sys.argv[1])
    out_file = Path(sys.argv[2])
    merge_fasta_files(src_dir, out_file)
