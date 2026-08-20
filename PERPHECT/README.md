# PERPHECT — PBI-Scope Integration

PERPHECT is a dual-CNN phage-host interaction predictor. This integration replaces PERPHECT's original CSV-based data loading with PBI-Scope's streaming database, enabling training on the full dataset without manual exports.

## What PBI-Scope Provides

PBI-Scope (`src/pbi/`) is a self-contained phage-host data platform built on DuckDB + pyfaidx. For ML training, it provides:

| PBI-Scope Component | What It Does |
|---------------------|--------------|
| `SequenceRetriever` | Core data engine. Connects DuckDB metadata to indexed FASTA files. Fetches phage/host sequences on demand via `get_phage_sequence()` and `get_host_sequence()`. Handles private data routing automatically. |
| `NegativeExampleGenerator` | Generates synthetic negative (non-interacting) pairs using random, GC-content, and taxonomy-based strategies. |
| `private_data` | Validates and ingests user-supplied phage-host interaction datasets. Manages `private_interactions` table with interaction type labels. |

### What We Did NOT Have to Build

PBI-Scope handled these natively — no custom code needed:

- DuckDB connection and query engine
- FASTA indexing and lazy sequence retrieval (pyfaidx, LRU cache, background preload)
- Private phage/host data routing and ingestion
- Interaction type storage and querying
- Negative example generation strategies
- Multi-contig genome assembly

## What the Adapter Handles

The `PBIAdapter` class (`pbi_adapter.py`) is the translation layer between PBI-Scope's data format and PERPHECT's expectations. **Everything here is custom code** — PBI-Scope provides no built-in PERPHECT support.

| Translation Layer | Why It's Needed |
|-------------------|-----------------|
| **String → Integer ID mapping** | PBI-Scope uses string IDs (`"GCF_000005845"`). PERPHECT's embedding layer requires contiguous integer indices. |
| **Length filtering** | PERPHECT expects bacteria ≤7M bp and phages ≤200K bp. Sequences outside these ranges are dropped. |
| **Zero-padding + truncation** | PERPHECT requires fixed-length inputs. Sequences are padded with zeros or truncated to the threshold. |
| **One-hot encoding** | PERPHECT expects `(length, 4)` numpy arrays (A/T/G/C channels). PBI-Scope returns raw strings. |
| **DataFrame schema mapping** | PERPHECT's original CSV format uses three tables (couples, bacteria, phages). The adapter produces these from PBI-Scope's pair DataFrames. |
| **TF generator interface** | PERPHECT's `model.fit()` expects an infinite-loop generator yielding `([bact_batch, phage_batch], labels)`. PBI-Scope has no Keras generator. |
| **Two-phase loading** | Loading all sequences up front takes 10+ minutes. The adapter queries pair IDs first (fast SQL), classifies them, applies limits, then fetches sequences only for selected pairs. |
| **Interaction classification** | The `private_interactions` table stores interaction types as strings. The adapter classifies `"no interaction"` / `"none"` / `"negative"` as label=0, everything else as label=1. |

### Adapter Usage

```python
from pbi_adapter import PBIAdapter

adapter = PBIAdapter(retriever, bacterium_threshold=7_000_000, phage_threshold=200_000)

# Phase 1: Query pair IDs only (fast SQL, no sequences loaded)
all_pairs = adapter.get_pair_ids_only(shuffle=True)

# Phase 2: Classify by interaction type
positive_pairs, negative_pairs = adapter.classify_pairs_by_interaction(all_pairs)

# Phase 3: Fetch sequences only for selected pairs (lazy)
couples, labels, sources = adapter.prepare_training_data(positive_pairs, negative_pairs)

# Phase 4: Create TF generator
gen = adapter.create_tf_generator(couples, labels, batch_size=4)
```

## Quick Start

### Prerequisites

1. PBI-Scope pipeline must have been run to build the database and sequence files
2. Docker with NVIDIA Container Toolkit (for GPU training)

### 1. Prepare Test Set (notebook)

Run once to create a held-out test set that is never used during training:

```bash
docker compose up -d analysis
# Open http://localhost:8886, navigate to PERPHECT/, open 01_prepare_test_set.ipynb
```

Outputs `test_data/test_set.csv` and `test_data/test_set.npz`.

### 2. Train (script)

```bash
# Quick test (1000 pairs, 3 epochs, excluding held-out test pairs)
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --limit 1000 --epochs 3 \
    --exclude-ids /workspace/PERPHECT/test_data/excluded_pairs.csv

# Full training
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --config /workspace/PERPHECT/config.yaml \
    --exclude-ids /workspace/PERPHECT/test_data/excluded_pairs.csv

# With cross-validation (5 folds)
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --limit 1000 --epochs 3 \
    --cross-validate 5 \
    --exclude-ids /workspace/PERPHECT/test_data/excluded_pairs.csv
```

### 3. Evaluate (notebook)

Open `02_evaluate_model.ipynb` to load the trained model and test set, then display:
- Classification report (precision, recall, F1)
- Confusion matrix
- ROC-AUC curve
- Precision-recall curve
- Prediction distribution histogram
- Per-source metrics (private_data vs generated negatives)

Results are saved to `./outputs/` on the host.

## Training Script Reference

### Command-Line Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--epochs` | 10 | Training epochs |
| `--batch-size` | 4 | Batch size |
| `--steps-per-epoch` | 400 | Batches per epoch |
| `--patience` | 5 | Early stopping patience |
| `--learning-rate` | 0.0004 | Initial learning rate |
| `--limit` | None | Limit positive pairs (None = all) |
| `--negative-ratio` | 1.0 | Negatives per positive |
| `--bacterium-threshold` | 7000000 | Max bacteria sequence length |
| `--phage-threshold` | 200000 | Max phage sequence length |
| `--bacterium-min-length` | 150000 | Min bacteria length to keep |
| `--phage-min-length` | 1500 | Min phage length to keep |
| `--output-dir` | /results | Output directory (mapped to `./outputs` on host) |
| `--run-name` | timestamp | Run name |
| `--config` | None | YAML config file |
| `--cross-validate` | 0 | K folds for stratified K-fold CV (0 = standard split) |
| `--exclude-ids` | None | CSV with Phage_ID,Host_ID to exclude from training (prevents data leakage) |
| `--no-gpu` | False | Force CPU |
| `--gpu-device` | 0 | GPU device index |
| `--verbose` | False | Verbose logging |
| `--log-file` | None | Log to file |

### Configuration File

```yaml
training:
  epochs: 20
  batch_size: 4
  steps_per_epoch: 800
  patience: 10
  learning_rate: 0.0004

data:
  negative_ratio: 1.0
  bacterium_threshold: 7000000
  phage_threshold: 200000

output:
  dir: /results

gpu:
  enabled: true
```

### Output Files

| File | Description |
|------|-------------|
| `model_best.keras` | Best model (by validation loss) |
| `model_final.keras` | Final model at end of training |
| `training_log.csv` | Epoch-by-epoch metrics |
| `accuracy.png` | Accuracy plot |
| `val_loss.png` | Loss plot |
| `results_test_set.csv` | Test set predictions |
| `summary.json` | Run metadata and metrics |
| `config.json` | Parameters used for this run |

## Architecture

```
Bacteria Input (7M bp, 4 channels)
    ↓
Conv1D(64, k=30, s=10) → MaxPool(15, s=5)
    ↓
Conv1D(32, k=25, s=10) → MaxPool(10, s=5)
    ↓
Conv1D(32, k=10, s=5) → MaxPool(2, s=2)
    ↓
Flatten → bacteria_features

Phage Input (200K bp, 4 channels)
    ↓
Conv1D(64, k=30, s=10) → MaxPool(15, s=5)
    ↓
Conv1D(32, k=25, s=10) → MaxPool(2, s=2)
    ↓
Flatten → phage_features

[bacteria_features | phage_features] → Dense(100) → Dropout(0.1) → Dense(1, sigmoid)
```

## Data Pipeline

### Pair Retrieval

1. **Shuffled ID query** (`get_pair_ids_only(shuffle=True)`): Fast SQL query using `ORDER BY MD5(Phage_ID || Host_ID)` for deterministic shuffling without loading sequences.
2. **Interaction classification** (`classify_pairs_by_interaction()`): Queries `private_interactions` to label each pair. Pairs with `"no interaction"` / `"none"` / `"negative"` become negatives (label=0). All others are positives (label=1).
3. **Limit applied to positives only**: Private-data negatives are always included.

### Negative Handling

| Source | Description | Label |
|--------|-------------|-------|
| `positive` | Known interacting pair from `phage_host_associations` | 1 |
| `private_data` | True negative from `private_interactions` (interaction = "no interaction") | 0 |
| `generated` | Synthetic negative from `NegativeExampleGenerator` | 0 |

Generated negatives are:
- **Deduplicated** against private-data negatives (no duplicate pairs)
- **Capped** at 2x the count of private-data negatives (prevents synthetic data from overwhelming true negatives)

### Train/Test/Validation Split

The data is split with **3-group stratification** to ensure each split contains a proportional mix of positives, private-data negatives, and generated negatives:

```
stratify_key = "pos" | "neg_private_data" | "neg_generated"
```

This ensures:
- The positive/negative ratio is preserved across splits
- True negatives from private data appear in all splits (not concentrated in one)
- Generated negatives are evenly distributed

Split sizes: 70% train / 15% validation / 15% test (stratified, shuffled, `random_state=42`).

### Split Source Distribution (logged during training)

```
train: {'positive': N, 'private_data': M, 'generated': K}
valid: {'positive': N', 'private_data': M', 'generated': K'}
test:  {'positive': N'', 'private_data': M'', 'generated': K''}
```

### Summary JSON

The `summary.json` file includes `split_sources` with per-split counts for each negative source, enabling analysis of model performance by data origin.

## Cross-Validation

Use `--cross-validate K` to run stratified K-fold cross-validation instead of a single train/val/test split:

```bash
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --cross-validate 5 --epochs 10
```

Each fold:
- Gets its own `fold_<k>/` directory with `model_best.keras`, `training_log.csv`, `summary.json`
- Uses 80/20 train/val split within the fold (stratified)
- Is evaluated on the held-out fold test set

After all folds, an aggregate `cv_summary.json` is saved with mean +/- std for val_loss, val_accuracy, and epochs trained.

Use the held-out test set from `01_prepare_test_set.ipynb` for final evaluation with `02_evaluate_model.ipynb`.

## Files

| File | Description |
|------|-------------|
| `train.py` | Production training script (CLI, GPU, YAML config, K-fold CV) |
| `pbi_adapter.py` | Adapter bridging PBI-Scope to PERPHECT format |
| `transforms.py` | One-hot encoding utilities |
| `plotting_utils.py` | Training visualization utilities |
| `config.yaml` | Default training configuration |
| `01_prepare_test_set.ipynb` | Create held-out test set (CSV + NPZ) |
| `02_evaluate_model.ipynb` | Load model, predict test set, display metrics |

## Testing

```bash
python -m pytest tests/test_perpherct_*.py -v
```

## Troubleshooting

### cuDNN on Pascal GPUs (GTX 1080 Ti)

cuDNN 9.x (bundled with TF 2.16+) dropped Pascal (sm_6.1) support. The Docker image pins TensorFlow to 2.15.x with cuDNN 8.9.x. After changing the Dockerfile, rebuild:

```bash
docker compose build --no-cache analysis
```

### Out of memory

- Reduce `--batch-size` (try 2)
- Reduce `--bacterium-threshold` and `--phage-threshold`
- Use `--limit` to train on a subset first

### Slow training

- Verify GPU is detected in logs (`GPU detected: 1 device(s)`)
- Check `nvidia-smi` on the host
- Reduce sequence thresholds for faster iteration

## Suggestions for PBI-Scope Adaptation

If you want to integrate another ML model with PBI-Scope, here are the key patterns from this integration:

### 1. Use the Two-Phase Loading Pattern

Don't load all sequences up front. PBI-Scope's `get_phage_host_pairs()` fetches everything — but for training you typically need a subset. Query IDs first, classify/filter, then fetch sequences only for selected pairs.

```python
# Instead of this (slow — fetches all sequences):
pairs = retriever.get_phage_host_pairs()

# Do this (fast — IDs only, then lazy fetch):
adapter = PBIAdapter(retriever, ...)
all_pairs = adapter.get_pair_ids_only()
positive, negative = adapter.classify_pairs_by_interaction(all_pairs)
couples, labels = adapter.prepare_training_data(positive, negative)
```

### 2. Handle ID Translation

PBI-Scope uses biological string IDs. Most ML frameworks expect integer indices. Build a bidirectional mapping:

```python
id_map = {string_id: int_idx for int_idx, string_id in enumerate(unique_ids)}
reverse_map = {v: k for k, v in id_map.items()}
```

### 3. Use the Negative Example Generator

PBI-Scope's `NegativeExampleGenerator` provides three strategies (random, GC-based, taxonomy-based). For interaction prediction, the `"mixed"` strategy is recommended. You can also use the `private_interactions` table directly — pairs labeled `"no interaction"` are true negatives.

### 4. Lazy Sequence Retrieval

PBI-Scope caches sequences in memory (LRU, up to 1000 entries). For training with generators, fetch sequences on demand rather than pre-loading. This keeps memory usage bounded by batch size, not dataset size.

### 5. Custom CNN Architectures

PBI-Scope does not enforce any model architecture. The adapter handles data formatting — you can plug in any model that accepts `(length, 4)` one-hot encoded inputs. For protein-level models, use `get_protein_sequences()` instead and adjust the encoding dimension.

### 6. Private Data Integration

Private datasets ingested via `private_data.ingest_private_sources_into_db()` are automatically visible to `SequenceRetriever`. The adapter's `classify_pairs_by_interaction()` queries the `private_interactions` table to distinguish positive/negative pairs. No extra work needed.

### 7. PyTorch Alternative

PBI-Scope provides PyTorch `Dataset` and `IterableDataset` classes (`create_streaming_dataset()`, `create_indexed_dataset()`). If your model uses PyTorch instead of Keras, use these directly — they handle batching, transforms, and private data routing internally.
