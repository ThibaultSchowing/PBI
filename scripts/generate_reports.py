#!.pixi/envs/reporting/bin/python
# Generate HTML reports from the merged CSV files 

import sys
import os
import pandas as pd
import random
import logging
from ydata_profiling import ProfileReport
import psutil

logging.basicConfig(level=logging.INFO)

def reservoir_sampling(file_path, sample_size):
    """
    Reservoir sampling: sample sample_size rows from a file sequentially and uniformly. 
    This function reads the file in chunks to handle large files efficiently and 
    returns a DataFrame containing the sampled rows.
    Too slow but can be used for smaller files.
    """
    sample = []
    reader = pd.read_csv(file_path, chunksize=1)
    logging.info(f"Sampling {sample_size} rows from {file_path} using reservoir sampling....")
    for i, row in enumerate(reader, start=1):
        if i <= sample_size:
            sample.append(row)
        else:
            j = random.randint(1, i)
            if j <= sample_size:
                sample[j-1] = row
    logging.info(f"Sampled {len(sample)} rows from the file. Now concatenating them into a DataFrame.")
    df_sample = pd.concat(sample, ignore_index=True)
    return df_sample

def fast_sample_known_size(file_path, sample_size, total_rows):
    import numpy as np
    logging.info(f"Sampling {sample_size} rows from {file_path} using fast sampling with known size....")
    sample_indices = set(random.sample(range(total_rows), sample_size))
    sample = []
    current_row = 0

    chunksize = 100000  # 100k lignes par chunk

    for chunk in pd.read_csv(file_path, chunksize=chunksize):
        chunk_indices = range(current_row, current_row + len(chunk))
        rows_to_keep = list(sample_indices.intersection(chunk_indices))
        if rows_to_keep:
            sample_rows = [i - current_row for i in rows_to_keep]
            sample.append(chunk.iloc[sample_rows])
        current_row += len(chunk)
        if len(set().union(*[df.index for df in sample])) >= sample_size:
            break

    df_sample = pd.concat(sample, ignore_index=True)
    return df_sample

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    logging.info(f"Input: {input_file}")
    logging.info(f"Output: {output_file}")

    # Check if input file exists
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} does not exist.")
        sys.exit(1)

    # Check file size in bytes
    file_size_bytes = os.path.getsize(input_file)
    file_size_mb = file_size_bytes / (1024 ** 2)
    logging.info(f"File size: {file_size_mb:.2f} MB")

    sample_size = 50000  # Adjust as needed

    try:
        # SMALL FILE: Read the entire file and sample
        if file_size_mb < 500:
            # Small file: read entirely and sample
            logging.info("File is smaller than 500 MB, reading entirely.")
            df = pd.read_csv(input_file)
            if df.empty:
                raise ValueError("The input DataFrame is empty. Cannot generate report.")
            sample_size_adjusted = min(len(df), sample_size)
            df_sampled = df.sample(n=sample_size_adjusted, random_state=42)
        # LARGE FILE: Use memory size checks to select sampling method
        else:
            # Large file, large memory: if memory is higher than 40 GB, read the entire file. Still faster than sampling method.
            if psutil.virtual_memory().available > 40 * 1024**3:
                logging.info("Available memory is sufficient, reading the entire file.")
                df = pd.read_csv(input_file)
                if df.empty:
                    raise ValueError("The input DataFrame is empty. Cannot generate report.")
                sample_size_adjusted = min(len(df), sample_size)
                df_sampled = df.sample(n=sample_size_adjusted, random_state=42)
            else:
            # Large file, low memory: use alternative sampling method
                logging.info("File is larger than 500 MB, using alternative sampling.")
                # Too slow !
                #df_sampled = reservoir_sampling(input_file, sample_size)

                # We know that the data has 43088582 rows, so we can use a faster method # !!! maybe not all files have the same number of rows
                #total_rows = 43088582 - 1  # Known total rows in the dataset
                # Get the total rows from the file without loading the entire file in memory (use sed or wc)
                total_rows = sum(1 for _ in open(input_file)) - 1
                logging.info(f"Total rows in the file: {total_rows}")
                df_sampled = fast_sample_known_size(input_file, sample_size, total_rows)

            if df_sampled.empty:
                raise ValueError("The input DataFrame is empty after sampling. Cannot generate report.")

        # Generate profiling report
        profile = ProfileReport(df_sampled, title=f"Profiling Report of {os.path.basename(input_file)}", explorative=True)
        profile.to_file(output_file=output_file)
        logging.info("Report generation completed successfully.")

    except Exception as e:
        logging.error(f"Error during report generation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
