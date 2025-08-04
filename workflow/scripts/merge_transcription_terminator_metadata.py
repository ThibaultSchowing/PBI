#!.pixi/envs/default/bin/python

import sys
print(f"Using python from: {sys.executable}")
import pandas as pd
import os
import logging
logging.basicConfig(level=logging.INFO)

# Snakemake inputs and outputs
inputs = snakemake.input
output = snakemake.output[0]

# Liste des DataFrames
dfs = []

for infile in inputs:
    print(f"Processing file: {infile}")
    df = pd.read_csv(infile, sep="\t")
    
    # Ajoute la colonne 'Phage_source' si elle n'existe pas
    if 'Phage_source' not in df.columns:
        source_name = os.path.basename(infile).split("_")[0]
        print(f"Source name derived from file: {source_name}")
        # Log info the source_name
        logging.info(f"Processing {infile} with source name '{source_name}'")
        df['Phage_source'] = source_name
        logging.info(f"Added 'Phage_source' column to {infile} with value '{source_name}'")

    dfs.append(df)
    

# Concatène tous les DataFrames
logging.info(f"Merging {len(dfs)} DataFrames")
merged_df = pd.concat(dfs, ignore_index=True)
logging.info(f"Merged DataFrame shape: {merged_df.shape}")

# Crée le dossier output si besoin
logging.info(f"Creating output directory if it does not exist: {os.path.dirname(output)}")
os.makedirs(os.path.dirname(output), exist_ok=True)

# Sauvegarde
logging.info(f"Saving merged DataFrame to {output}")
merged_df.to_csv(output, index=False)

print(f"[INFO] Merged {len(inputs)} files into {output} with shape {merged_df.shape}")
