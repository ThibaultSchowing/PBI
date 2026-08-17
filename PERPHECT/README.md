# PERPHECT

PERPHECT is a phage-host (phage-bacteria) interaction predictor based on a dual CNN architecture. It processes zero-padded, one-hot-encoded DNA sequences of bacteria and phages to predict whether a given phage can infect a given bacterium.

This integration replaces PERPHECT's original CSV-based data loading with PBI-Scope's streaming data API, enabling training on the full database without pre-computing fixed-length padded CSV files.

## Architecture

PERPHECT uses two separate CNN branches:
- **Bacteria branch**: 3 convolutional blocks processing 7M bp sequences
- **Phage branch**: 2 convolutional blocks processing 200K bp sequences
- **Classification head**: Concatenated features → Dense(100) → Dropout(0.1) → Dense(1, sigmoid)

## Data Format

The model expects data in these formats:

| DataFrame | Columns | Description |
|-----------|---------|-------------|
| Phages | `phage_id, phage_sequence` | Phage sequences (one-hot encoded) |
| Bacteria | `bacterium_id, bacterium_sequence` | Bacteria sequences (one-hot encoded) |
| Couples | `id, bacterium_id, phage_id, interaction_type` | Interaction pairs (0 = negative, 1 = positive) |

### Sequence Thresholds

| Parameter | Value | Description |
|-----------|-------|-------------|
| `BACTERIUM_THRESHOLD` | 7,000,000 bp | Maximum bacteria sequence length (zero-padded) |
| `PHAGE_THRESHOLD` | 200,000 bp | Maximum phage sequence length (zero-padded) |
| `BACTERIUM_MIN_LENGTH` | 150,000 bp | Minimum bacteria length to include |
| `PHAGE_MIN_LENGTH` | 1,500 bp | Minimum phage length to include |

## Quick Start

### Prerequisites

1. The PBI-Scope pipeline must have been run to build the database
2. Docker must be installed with the `docker-compose.yml` from the project root

### 1. Start the Analysis Container

From the project root directory:

```bash
docker compose up -d analysis
```

### 2. Connect to JupyterLab

Create an SSH tunnel to the analysis container:

```bash
ssh -L 8886:localhost:8888 <your-host>
```

Then open `http://localhost:8888` in your browser (or `http://localhost:8886` if using SSH tunnel).

If you configured a token in `docker-compose.yml` (the `JUPYTER_TOKEN` environment variable), enter it when prompted.

### 3. Navigate to PERPHECT

In JupyterLab's file browser, navigate to the `PERPHECT/` folder.

### 4. Open the Training Notebook

Open `PBI_Perphect_Training.ipynb` and follow the instructions.

## Files

| File | Description |
|------|-------------|
| `PBI_Perphect_Training.ipynb` | Main training notebook (start here) |
| `model_with_different_paddings_mixed.py` | Original PERPHECT model (standalone) |
| `pbi_adapter.py` | Adapter bridging PBI-Scope to PERPHECT format |
| `transforms.py` | One-hot encoding utilities |
| `plotting_utils.py` | Training visualization utilities |
| `Notebook.ipynb` | Original exploratory notebook |

## How the Adapter Works

The `PBIAdapter` class (`pbi_adapter.py`) handles the translation between PBI-Scope's data format and PERPHECT's expectations:

1. **ID Mapping**: PBI-Scope uses string IDs (e.g., `GCF_000012345.1`), while PERPHECT uses integer IDs. The adapter maintains bidirectional mappings.

2. **Sequence Fetching**: Sequences are loaded on-demand from PBI-Scope's DuckDB database and indexed FASTA files.

3. **Length Filtering**: Sequences shorter than the minimum thresholds are excluded.

4. **Padding**: Sequences are zero-padded to fixed lengths (or truncated if longer).

5. **One-Hot Encoding**: DNA sequences are converted to numpy arrays of shape `(length, 4)` with columns `[A, T, G, C]`.

### DataFrame Mode

For small datasets or validation:

```python
from pbi_adapter import PBIAdapter

adapter = PBIAdapter(retriever)
couples_df, bacteria_df, phages_df = adapter.to_perphect_dataframes(
    positive_pairs, negative_pairs
)
```

### Generator Mode

For memory-efficient training on large datasets:

```python
# Prepare training arrays
couples, labels = adapter.prepare_training_data(positive_pairs, negative_pairs)

# Create TF-compatible generator
gen = adapter.create_tf_generator(couples, labels, batch_size=16)

# Use with model.fit()
model.fit(gen, steps_per_epoch=400, epochs=10)
```

## Training on Full Data

To train on all available data instead of a subset:

```python
# Remove the limit parameter
positive_pairs = retriever.get_phage_host_pairs(host_contig_mode='concat')

# Generate negatives
negative_pairs = neg_gen.generate_random_negatives(positive_pairs, ratio=1.0)

# Prepare and train
couples, labels = adapter.prepare_training_data(positive_pairs, negative_pairs)
```

**Note:** Training on the full database requires sufficient RAM for sequence caching. The adapter caches sequences in memory as they are fetched.

## Running Tests

```bash
# Run PERPHECT-specific tests
python -m pytest tests/test_perpherct_transforms.py tests/test_perpherct_plotting.py tests/test_perpherct_adapter.py -v

# Run all tests
python -m pytest tests/ -v
```

## Original PERPHECT (Standalone)

The original PERPHECT model can be run independently with pre-computed CSV files:

```bash
python model_with_different_paddings_mixed.py --batch 16 --epoch 3
```

This requires CSV files in `data/mixed_data_set_/` with the format described in the Data Format section above.
