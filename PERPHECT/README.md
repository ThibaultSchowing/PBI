# PERPHECT — Phage-Host Interaction Predictor

PERPHECT is a dual-CNN model that predicts phage-bacteria interactions from genomic sequences. It uses two parallel convolutional encoders — one for the bacterial genome (up to 7M bp) and one for the phage genome (up to 200K bp) — whose features are concatenated and fed to a classification head.

This directory contains the complete training pipeline, integrated with PBI-Scope as the data backend.

## How It Works

### Training Pipeline Overview

```
PBI-Scope DB (DuckDB + FASTA)
    │
    ▼
PBIAdapter.get_pair_ids_only()     ← fast SQL, no sequences
    │
    ▼
PBIAdapter.classify_pairs_by_interaction()
    │
    ▼
NegativeExampleGenerator           ← synthetic negatives
    │
    ▼
PBIAdapter.prepare_training_data() ← fetches sequences, one-hot encodes
    │
    ▼
train_test_split()                 ← 70/15/15 stratified
    │
    ▼
model.fit()                        ← focal loss, val_auc monitoring
    │
    ▼
model_best.keras + summary.json
```

### The Adapter (`pbi_adapter.py`)

The `PBIAdapter` translates between PBI-Scope's data format and PERPHECT's expectations:

| What It Does | Why |
|---|---|
| **String → integer ID mapping** | PBI-Scope uses `"GCF_000005845"`. PERPHECT's embedding layer needs contiguous integers. |
| **Length filtering** | Drops bacteria >7M bp or <150K bp. Drops phages >200K bp or <1.5K bp. |
| **Zero-padding + truncation** | PERPHECT requires fixed-length `(length, 4)` inputs. |
| **One-hot encoding** | Converts `"ATGC..."` strings to `(length, 4)` numpy arrays via a 256-entry LUT. |
| **LRU host cache** | `OrderedDict` with max 200 entries (~5.6GB cap). Hosts are reused across epochs; phages are re-encoded on the fly (~5ms). |
| **Interaction classification** | Queries `private_interactions` table. Labels: `"no interaction"` / `"none"` / `"negative"` → 0, everything else → 1. |
| **TF generator** | Yields `([bact_batch, phage_batch], labels)` indefinitely for `model.fit()`. |

### Data Sources

| Source | Label | Description |
|---|---|---|
| `positive` | 1 | Known interacting pair from `phage_host_associations` |
| `private_data` | 0 | True negative from `private_interactions` |
| `generated` | 0 | Synthetic negative from `NegativeExampleGenerator` |

Generated negatives are deduplicated against private-data negatives and capped at `--negative-ratio` × private negative count.

### Stratification

Data is split with a 3-group stratification key to ensure proportional distribution:

```
stratify_key = "pos" | "neg_private_data" | "neg_generated"
```

Default split: 70% train / 15% validation / 15% test.

## Quick Start

### Prerequisites

1. PBI-Scope pipeline must have been run to build the database and sequence files
2. Docker with NVIDIA Container Toolkit (for GPU training)

### 1. Prepare Test Set

Run once to create a held-out test set that is never used during training:

```bash
docker compose up -d analysis
# Open http://localhost:8886 → PERPHECT/01_prepare_test_set.ipynb
```

Outputs: `test_data/test_set.csv`, `test_data/test_set.npz`, `test_data/excluded_pairs.csv`

### 2. Train

```bash
# Quick test (1000 pairs, 3 epochs)
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --limit 1000 --epochs 3

# Pre-training (excludes private data, 5-fold CV)
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --config /workspace/PERPHECT/config.yaml \
    --profile pretrain --gpu-device 3

# Fine-tuning (requires pre-trained model from stage 1)
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --config /workspace/PERPHECT/config.yaml \
    --profile finetune \
    --pretrained-model /workspace/PERPHECT/outputs/run_*/fold_1/model_best.keras \
    --freeze-base --gpu-device 3
```

### 3. Evaluate

Open `02_evaluate_model.ipynb` to load a trained model against the held-out test set.

## Two-Stage Workflow

For best results with private data:

| Stage | Profile | Data Used | Output |
|---|---|---|---|
| Pre-train | `--profile pretrain` | All PBI-Scope **except** private source | `model_best.keras` |
| Fine-tune | `--profile finetune` | **Only** private source | `model_finetuned_best.keras` |

Key rules:
- `--exclude-ids` (the held-out test set) is excluded in **both** stages
- Fine-tuning uses `--freeze-base` to freeze CNN encoders, training only the classification head
- Fine-tuning uses lower LR (0.0001), fewer epochs (5), smaller batch (16)

## Model Architecture

```
Bacteria Input (7,000,000 bp × 4 channels)
    │
    Conv1D(64, k=30, s=10) → MaxPool(15, s=5)
    Conv1D(32, k=25, s=10) → MaxPool(10, s=5)
    Conv1D(32, k=10, s=5)  → MaxPool(2, s=2)
    │
    Flatten → bacteria_features

Phage Input (200,000 bp × 4 channels)
    │
    Conv1D(64, k=30, s=10) → MaxPool(15, s=5)
    Conv1D(32, k=25, s=10) → MaxPool(2, s=2)
    │
    Flatten → phage_features

Concatenate → Dense(100, relu) → Dropout(0.1) → Dense(1, sigmoid)
```

Bacteria has 3 conv layers (larger input), phage has 2 conv layers (smaller input).

## Configuration

### Config File (`config.yaml`)

Uses a `defaults` + `profiles` structure. Priority: **CLI args > profile overrides > defaults**.

```yaml
defaults:
  training:
    epochs: 15
    batch_size: 32
    patience: 5
    learning_rate: 0.0004
    focal_alpha: 0.25
    focal_gamma: 2.0
    cross_validate: 0

  data:
    negative_ratio: 1.0
    bacterium_threshold: 7000000
    phage_threshold: 200000
    bacterium_min_length: 150000
    phage_min_length: 1500
    exclude_ids: "test_data/excluded_pairs.csv"
    exclude_sources:
      - "PERPHECT_private"

  output:
    dir: /results

  gpu:
    enabled: true

profiles:
  pretrain:
    training:
      cross_validate: 5

  finetune:
    training:
      epochs: 5
      batch_size: 16
      patience: 3
      fine_tune_lr: 0.0001
      fine_tune_epochs: 5
      freeze_base: true
      freeze_up_to: "concatenated_features"
      finetuned_model_name: "model_finetuned_best.keras"
```

Parameters like `exclude_ids` and `exclude_sources` are config-only — not CLI flags.

### Command-Line Arguments

| Argument | Default | Description |
|---|---|---|
| `--epochs` | 15 | Training epochs |
| `--batch-size` | 32 | Batch size |
| `--steps-per-epoch` | None | Batches per epoch (None = full dataset) |
| `--patience` | 5 | Early stopping patience |
| `--learning-rate` | 0.0004 | Initial learning rate |
| `--limit` | None | Limit positive pairs (None = all) |
| `--negative-ratio` | 1.0 | Max synthetic negatives × private negatives |
| `--focal-alpha` | 0.25 | Focal loss positive-class weight |
| `--focal-gamma` | 2.0 | Focal loss focusing parameter |
| `--bacterium-threshold` | 7000000 | Max bacteria length |
| `--phage-threshold` | 200000 | Max phage length |
| `--bacterium-min-length` | 150000 | Min bacteria length |
| `--phage-min-length` | 1500 | Min phage length |
| `--output-dir` | /results | Output directory (Docker: mapped to `./outputs`) |
| `--run-name` | timestamp | Run name |
| `--config` | None | YAML config file |
| `--profile` | None | Profile name (`pretrain`, `finetune`) |
| `--cross-validate` | 0 | K folds for CV (0 = standard split) |
| `--pretrained-model` | None | Path to pre-trained `.keras` model |
| `--freeze-base` | False | Freeze CNN encoders |
| `--freeze-up-to` | concatenated_features | Layer to freeze up to |
| `--fine-tune-lr` | 0.0001 | Fine-tuning learning rate |
| `--fine-tune-epochs` | 5 | Fine-tuning epochs |
| `--finetuned-model-name` | model_finetuned_best.keras | Fine-tuned model output name |
| `--no-gpu` | False | Force CPU |
| `--gpu-device` | 0 | GPU device index |

## Loss Function & Metrics

### Focal Loss

Training uses **focal loss** (Lin et al., 2017) instead of binary cross-entropy. Focal loss down-weights easy negatives automatically, so the model focuses on hard examples. This is critical because generated negatives (random phage-host pairs) are trivially distinguishable from real interactions — standard BCE lets the model achieve high accuracy by memorizing easy negatives.

Parameters: `--focal-alpha` (class balance, default 0.25) and `--focal-gamma` (focusing, default 2.0).

### Metrics

| Metric | When | Purpose |
|---|---|---|
| `val_auc` | Training | Early stopping + checkpoint selection (threshold-independent) |
| `accuracy` | Training | Logged for reference |
| MCC | Test | Matthews Correlation Coefficient — balanced even with class imbalance |
| Sensitivity | Test | Recall for positive class |
| Specificity | Test | Recall for negative class |
| Classification report | Test | Precision, recall, F1 per class |

## Output Files

### Standard Split

```
outputs/run_<timestamp>/
├── model_best.keras              # Best model (by val_auc)
├── model_final.keras             # Final model
├── training_log.csv              # Epoch metrics
├── accuracy.png                  # Training/validation accuracy
├── val_loss.png                  # Training/validation loss
├── pairs_all.csv                 # All pairs (string IDs, labels, sources)
├── pairs_train.csv               # Training split pairs
├── pairs_val.csv                 # Validation split pairs
├── pairs_test.csv                # Test split pairs
├── results_test_set.csv          # Test predictions (string + integer IDs)
├── summary.json                  # Run metadata + metrics
└── config.json                   # Parameters used
```

### Cross-Validation

```
outputs/run_<timestamp>/
├── fold_1/
│   ├── model_best.keras
│   ├── training_log.csv
│   ├── pairs_train.csv
│   ├── pairs_val.csv
│   ├── pairs_test.csv
│   ├── summary.json
│   └── ...
├── fold_2/
│   └── ...
├── cv_summary.json               # Aggregated mean ± std across folds
└── ...
```

### Traceability

Every run saves `pairs_all.csv` and per-split `pairs_*.csv` with original string IDs (`host_id`, `phage_id`), labels, and data sources. The `results_test_set.csv` includes both string and integer IDs alongside predictions.

This lets you:
- Verify which specific pairs ended up in each split
- Check for data leakage between train and test
- Analyze per-species or per-source model performance
- Cross-reference with taxonomy metadata

## Cross-Validation

```bash
python train.py --config config.yaml --cross-validate 5 --epochs 10
```

- Each fold gets its own `fold_X/` directory with model, logs, and pair CSVs
- Uses 80/20 train/val split within each fold
- After all folds: `cv_summary.json` with mean ± std for val_loss, val_accuracy, and epochs
- **Not compatible** with `--pretrained-model` (fine-tuning)

## Evaluation Notebook

`02_evaluate_model.ipynb` loads a trained model and evaluates it against the held-out test set from `01_prepare_test_set.ipynb`:

1. Loads `test_set.npz` (pre-encoded numpy arrays)
2. Finds `model_best.keras` or `model_final.keras` in the outputs directory
3. Generates predictions
4. Computes: MCC, sensitivity, specificity, classification report
5. Plots: confusion matrix, ROC curve, precision-recall curve, prediction distribution
6. Per-source metrics (if sources are available)
7. Optionally compares multiple runs side by side

## Files

| File | Description |
|---|---|
| `train.py` | Training script (CLI, GPU, config profiles, K-fold CV, fine-tuning) |
| `pbi_adapter.py` | Adapter bridging PBI-Scope data to PERPHECT format |
| `transforms.py` | Vectorized one-hot encoding via numpy LUT |
| `plotting_utils.py` | Training history plots |
| `config.yaml` | Training configuration (defaults + profiles) |
| `01_prepare_test_set.ipynb` | Creates held-out test set (CSV + NPZ) |
| `02_evaluate_model.ipynb` | Loads model, evaluates on test set, compares runs |

## Testing

```bash
python -m pytest tests/test_perpherct_*.py
```

Test coverage:
- `test_perpherct_transforms.py` — one-hot encoding edge cases
- `test_perpherct_adapter.py` — ID mapping, encoding, generators, interaction classification
- `test_perpherct_model.py` — model building (skipped without keras)
- `test_perpherct_train.py` — CLI args, config loading, profile merging, full pretrain/finetune paths (mocked)
- `test_perpherct_plotting.py` — training history plots

## Troubleshooting

**cuDNN on Pascal GPUs (GTX 1080 Ti):** cuDNN 9.x dropped Pascal (sm_6.1) support. The Docker image pins TF 2.15 + cuDNN 8.9.x. Rebuild with `docker compose build --no-cache analysis`.

**Out of memory:** Reduce `--batch-size` (try 2), reduce thresholds, or use `--limit` on a subset.

**Slow training:** Check `nvidia-smi` on host. Verify GPU detected in logs (`GPU detected: 1 device(s)`). The host cache (~5.6GB) ensures each unique host is encoded once.
