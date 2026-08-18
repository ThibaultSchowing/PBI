# PERPHECT — PBI-Scope Integration Tutorial

PERPHECT is a phage-host interaction predictor using a dual CNN architecture. This integration demonstrates how to use PBI-Scope's data streaming to train PERPHECT on the full database, replacing the original CSV-based workflow.

This folder serves as a complete tutorial for adapting PBI-Scope data to your own ML models.

## Overview

| What | How |
|------|-----|
| **Explore interactively** | Jupyter notebook (`PBI_Perphect_Training.ipynb`) |
| **Train at scale** | CLI script (`train.py`) via `docker compose run` |
| **GPU acceleration** | NVIDIA CUDA container (automatic if GPU available) |

## Prerequisites

### 1. PBI-Scope Pipeline

The pipeline must have been run to build the database and sequence files:

```bash
docker compose up -d pipeline
# Wait for completion...
docker compose run --rm pipeline snakemake --cores 2 --snakefile /app/workflow/Snakefile
```

### 2. Docker with GPU Support (Optional but Recommended)

The analysis container includes NVIDIA CUDA for GPU-accelerated training. To use GPU:

**Install NVIDIA Container Toolkit on the host:**

```bash
# Ubuntu/Debian
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | \
  sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg

curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

**Verify GPU access:**

```bash
docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi
```

**Without GPU:** Training works on CPU but is significantly slower. The container runs fine without the NVIDIA toolkit — just skip the GPU verification.

## Quick Start

### Option 1: Jupyter Notebook (Exploration)

Best for: testing code, visualizing results, understanding the data pipeline.

```bash
# Start the analysis container
docker compose up -d analysis

# Connect via SSH tunnel
ssh -L 8886:localhost:8888 <your-host>

# Open http://localhost:8886, navigate to PERPHECT/, open the notebook
```

The notebook walks through:
1. Connecting to PBI-Scope
2. Loading positive pairs
3. Generating negative examples
4. Transforming data with the adapter
5. Building and training the model
6. Evaluating results

### Option 2: Training Script (Production)

Best for: long training runs, background execution, reproducible experiments.

```bash
# Quick test run (1000 pairs, 3 epochs)
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --limit 1000 --epochs 3

# Full training with config file
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --config /workspace/PERPHECT/config.yaml

# Custom parameters
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py \
    --epochs 20 --batch-size 32 --output-dir /results/my_run

# Force CPU even if GPU available
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --no-gpu
```

**Monitor training progress:**

```bash
# Follow container logs
docker compose logs -f analysis

# Check output directory
ls -la outputs/<run_name>/
```

### Option 3: Notebook-to-Script Conversion

For keeping the notebook and a runnable script in sync:

```bash
docker compose run --rm analysis \
  python /workspace/PERPHECT/notebook_to_script.py \
    --verify /workspace/PERPHECT/train.py
```

## Training Script Reference

### Command-Line Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--epochs` | 10 | Training epochs |
| `--batch-size` | 16 | Batch size |
| `--steps-per-epoch` | 400 | Batches per epoch |
| `--patience` | 5 | Early stopping patience |
| `--learning-rate` | 0.0004 | Initial learning rate |
| `--limit` | None | Limit positive pairs (None = all) |
| `--negative-ratio` | 1.0 | Negatives per positive |
| `--bacterium-threshold` | 7000000 | Max bacteria sequence length |
| `--phage-threshold` | 200000 | Max phage sequence length |
| `--bacterium-min-length` | 150000 | Min bacteria length to keep |
| `--phage-min-length` | 1500 | Min phage length to keep |
| `--output-dir` | results | Output directory |
| `--run-name` | timestamp | Run name |
| `--config` | None | YAML config file |
| `--no-gpu` | False | Force CPU |
| `--verbose` | False | Verbose logging |
| `--log-file` | None | Log to file |

### Configuration File

Use `config.yaml` to save parameter sets:

```yaml
training:
  epochs: 20
  batch_size: 32
  steps_per_epoch: 800
  patience: 10
  learning_rate: 0.0004

data:
  # limit: null  # null = all data
  negative_ratio: 1.0
  bacterium_threshold: 7000000
  phage_threshold: 200000

output:
  dir: results

gpu:
  enabled: true
```

```bash
docker compose run --rm analysis \
  python /workspace/PERPHECT/train.py --config /workspace/PERPHECT/config.yaml
```

### Output Files

Each training run creates a timestamped directory with:

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

PERPHECT uses two separate CNN branches:

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

## How the Adapter Works

The `PBIAdapter` class bridges PBI-Scope's data format to PERPHECT's expectations:

1. **ID Mapping**: PBI-Scope string IDs → PERPHECT integer IDs
2. **Sequence Fetching**: On-demand from DuckDB + FASTA files
3. **Length Filtering**: Excludes sequences below minimum thresholds
4. **Padding**: Zero-padding to fixed lengths (or truncation)
5. **One-Hot Encoding**: DNA → numpy array of shape `(length, 4)`

```python
from pbi_adapter import PBIAdapter

adapter = PBIAdapter(retriever, bacterium_threshold=7_000_000, phage_threshold=200_000)

# DataFrame mode (small data)
couples_df, bacteria_df, phages_df = adapter.to_perphect_dataframes(positives, negatives)

# Generator mode (large data)
couples, labels = adapter.prepare_training_data(positives, negatives)
gen = adapter.create_tf_generator(couples, labels, batch_size=16)
```

## Files

| File | Description |
|------|-------------|
| `PBI_Perphect_Training.ipynb` | Interactive training notebook |
| `train.py` | Standalone training script for production |
| `config.yaml` | Default training configuration |
| `pbi_adapter.py` | Adapter bridging PBI-Scope to PERPHECT format |
| `transforms.py` | One-hot encoding utilities |
| `plotting_utils.py` | Training visualization utilities |
| `notebook_to_script.py` | Notebook-to-script converter with verification |
| `model_with_different_paddings_mixed.py` | Original PERPHECT model (standalone) |
| `Notebook.ipynb` | Original exploratory notebook |

## Testing

```bash
# Run all PERPHECT tests (unit tests, no GPU/database needed)
python -m pytest tests/test_perpherct_*.py -v

# Run specific test suites
python -m pytest tests/test_perpherct_adapter.py -v    # Adapter tests
python -m pytest tests/test_perpherct_model.py -v      # Model build tests (requires Keras)
python -m pytest tests/test_perpherct_train.py -v      # Script tests
python -m pytest tests/test_perpherct_notebook_script.py -v  # Converter tests
```

## GPU vs CPU

| Aspect | CPU | GPU |
|--------|-----|-----|
| Training speed | ~1 epoch/hour | ~10-20x faster |
| Memory | Lower | Higher (GPU VRAM) |
| Setup | No extra setup | NVIDIA Container Toolkit required |
| Docker image | N/A (CUDA works on CPU) | ~3GB larger |

The container auto-detects GPU at startup. If no GPU is found, it falls back to CPU automatically.

## Troubleshooting

### "No GPU detected" but GPU is installed

1. Verify NVIDIA driver: `nvidia-smi`
2. Verify Container Toolkit: `docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi`
3. Restart Docker: `sudo systemctl restart docker`

### Training is very slow

- Check GPU detection in the notebook or script logs
- Reduce `--bacterium-threshold` and `--phage-threshold` for testing
- Use `--limit` to train on a subset first
- Increase `--batch-size` if GPU memory allows

### Out of memory

- Reduce `--batch-size`
- Reduce `--bacterium-threshold` (default 7M bp is large)
- The adapter caches sequences in memory; use `--limit` for large datasets

### Container fails to start

- Check Docker logs: `docker compose logs analysis`
- Verify the database exists: `ls data/processed/databases/`
- Check port availability: `lsof -i :8886`
