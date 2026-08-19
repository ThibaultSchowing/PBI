#!/usr/bin/env python3
"""
PERPHECT training script for PBI-Scope data.

Standalone script for training the PERPHECT phage-host interaction predictor
using data streamed from PBI-Scope. Designed for long-running training jobs
that don't require an open browser tab.

Usage:
    # Quick test run
    docker compose run --rm analysis python /workspace/PERPHECT/train.py \
        --limit 1000 --epochs 3

    # Full training with config file
    docker compose run --rm analysis python /workspace/PERPHECT/train.py \
        --config /workspace/PERPHECT/config.yaml

    # Override specific parameters
    docker compose run --rm analysis python /workspace/PERPHECT/train.py \
        --epochs 20 --batch-size 32 --output-dir /results/my_run

    # Force CPU even if GPU available
    docker compose run --rm analysis python /workspace/PERPHECT/train.py \
        --no-gpu --epochs 10
"""

import argparse
import json
import logging
import math
import os
import sys
from datetime import datetime
from pathlib import Path

# Disable cuDNN autotuning and XLA to prevent autotuning failures on some GPU configs
os.environ["TF_CUDNN_USE_AUTOTUNER"] = "0"
os.environ["CUDNN_USE_AUTOTUNER"] = "0"
os.environ["TF_XLA_FLAGS"] = "--tf_xla_enable_xla_devices=false"

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

def setup_logging(log_file=None, verbose=False):
    """Configure logging to console and optional file."""
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s - %(levelname)s - %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(level=level, format=fmt, datefmt=datefmt, handlers=handlers)


# ---------------------------------------------------------------------------
# GPU detection
# ---------------------------------------------------------------------------

def detect_gpu(gpu_device=0):
    """Detect and report GPU availability. Returns True if GPU is available."""
    try:
        import tensorflow as tf
        tf.config.optimizer.set_jit(False)
        gpus = tf.config.list_physical_devices("GPU")
        if gpus:
            logging.info(f"GPU detected: {len(gpus)} device(s)")
            for gpu in gpus:
                logging.info(f"  - {gpu.name}")
            if gpu_device < len(gpus):
                tf.config.set_visible_devices(gpus[gpu_device], "GPU")
                logging.info(f"Selected GPU {gpu_device}: {gpus[gpu_device].name}")
            else:
                logging.warning(f"GPU {gpu_device} not found, using GPU 0")
                tf.config.set_visible_devices(gpus[0], "GPU")
            return True
        else:
            logging.info("No GPU detected. Training will use CPU (slower).")
            return False
    except Exception as e:
        logging.warning(f"Could not detect GPU: {e}")
        return False


# ---------------------------------------------------------------------------
# Model builder
# ---------------------------------------------------------------------------

def build_model(bacterium_threshold=7_000_000, phage_threshold=200_000):
    """Build the PERPHECT dual-CNN model."""
    import keras
    from keras.models import Model
    from keras.layers import (
        Input, Conv1D, MaxPooling1D, Flatten, Dense, Dropout, Concatenate
    )

    # Bacteria branch
    input1 = Input(shape=(bacterium_threshold, 4), name="bacterial_input")
    conv1_1 = Conv1D(64, 30, strides=10, activation="relu", name="bacterial_conv_1")(input1)
    maxpool1_1 = MaxPooling1D(15, strides=5, name="bacterial_maxpool_1")(conv1_1)
    conv1_2 = Conv1D(32, 25, strides=10, activation="relu", name="bacterial_conv_2")(maxpool1_1)
    maxpool1_2 = MaxPooling1D(10, strides=5, name="bacterial_maxpool_2")(conv1_2)
    conv1_3 = Conv1D(32, 10, strides=5, activation="relu", name="bacterial_conv_3")(maxpool1_2)
    maxpool1_3 = MaxPooling1D(2, strides=2, name="bacterial_maxpool_3")(conv1_3)
    flatten_bact = Flatten(name="bacteria_features")(maxpool1_3)

    # Phage branch
    input2 = Input(shape=(phage_threshold, 4), name="phage_input")
    conv2_1 = Conv1D(64, 30, strides=10, activation="relu", name="phage_conv_1")(input2)
    maxpool2_1 = MaxPooling1D(15, strides=5, name="phage_maxpool_1")(conv2_1)
    conv2_2 = Conv1D(32, 25, strides=10, activation="relu", name="phage_conv_2")(maxpool2_1)
    maxpool2_2 = MaxPooling1D(2, strides=2, name="phage_maxpool_2")(conv2_2)
    flatten_phage = Flatten(name="phage_features")(maxpool2_2)

    # Classification head
    concat_features = Concatenate(name="concatenated_features")([flatten_bact, flatten_phage])
    dense1 = Dense(100, activation="relu")(concat_features)
    dropout1 = Dropout(0.10)(dense1)
    dense2 = Dense(1, activation="sigmoid")(dropout1)

    model = Model(name="Perphect", inputs=[input1, input2], outputs=dense2)
    return model


# ---------------------------------------------------------------------------
# Learning rate schedule
# ---------------------------------------------------------------------------

def step_decay(epoch, initial_lrate=0.0004, drop=0.5, epochs_drop=4.0):
    """Learning rate step decay schedule."""
    return initial_lrate * math.pow(drop, math.floor((1 + epoch) / epochs_drop))


# ---------------------------------------------------------------------------
# Main training function
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Train PERPHECT phage-host interaction predictor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Training parameters
    parser.add_argument("--epochs", type=int, default=10, help="Number of training epochs (default: 10)")
    parser.add_argument("--batch-size", type=int, default=4, help="Batch size (default: 4)")
    parser.add_argument("--steps-per-epoch", type=int, default=400, help="Batches per epoch (default: 400)")
    parser.add_argument("--patience", type=int, default=5, help="Early stopping patience (default: 5)")
    parser.add_argument("--learning-rate", type=float, default=0.0004, help="Initial learning rate (default: 0.0004)")

    # Data parameters
    parser.add_argument("--limit", type=int, default=None, help="Limit positive pairs (None = all)")
    parser.add_argument("--negative-ratio", type=float, default=1.0, help="Negatives per positive (default: 1.0)")
    parser.add_argument("--bacterium-threshold", type=int, default=7_000_000, help="Max bacteria sequence length")
    parser.add_argument("--phage-threshold", type=int, default=200_000, help="Max phage sequence length")
    parser.add_argument("--bacterium-min-length", type=int, default=150_000, help="Min bacteria length to keep")
    parser.add_argument("--phage-min-length", type=int, default=1_500, help="Min phage length to keep")

    # Output parameters
    parser.add_argument("--output-dir", type=str, default="results", help="Output directory (default: results)")
    parser.add_argument("--run-name", type=str, default=None, help="Run name (default: timestamp)")

    # Configuration
    parser.add_argument("--config", type=str, default=None, help="YAML config file (overrides CLI args)")
    parser.add_argument("--no-gpu", action="store_true", help="Force CPU even if GPU available")
    parser.add_argument("--no-cudnn", action="store_true", help="Disable cuDNN (needed for Pascal/sm_6.1 GPUs with cuDNN 9.x)")
    parser.add_argument("--gpu-device", type=int, default=0, help="GPU device index (default: 0, use 1 for second GPU)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument("--log-file", type=str, default=None, help="Log to file in addition to console")

    args = parser.parse_args()

    # Load config file if provided
    if args.config:
        import yaml
        with open(args.config) as f:
            config = yaml.safe_load(f)
        # Config values override defaults but not explicit CLI args
        for section in ["training", "data", "output", "gpu"]:
            if section in config:
                for key, value in config[section].items():
                    attr = key.replace("-", "_")
                    # Only apply if the CLI arg is at its default value
                    current = getattr(args, attr, None)
                    if current is not None:
                        setattr(args, attr, value)

    # Setup
    setup_logging(args.log_file, args.verbose)
    logging.info("=" * 60)
    logging.info("PERPHECT Training Script")
    logging.info("=" * 60)

    # GPU detection
    if args.no_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        logging.info("GPU disabled by user (--no-gpu)")
    else:
        if args.no_cudnn:
            os.environ["TF_USE_CUDNN"] = "0"
            logging.info("cuDNN disabled (--no-cudnn): using TF native convolutions")
        gpu_available = detect_gpu(args.gpu_device)

    # Generate run name
    run_name = args.run_name or datetime.now().strftime("run_%Y%m%d_%H%M%S")
    output_dir = Path(args.output_dir) / run_name
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Output directory: {output_dir}")

    # Save configuration
    config_path = output_dir / "config.json"
    with open(config_path, "w") as f:
        json.dump(vars(args), f, indent=2, default=str)
    logging.info(f"Configuration saved to {config_path}")

    # -----------------------------------------------------------------------
    # Data loading
    # -----------------------------------------------------------------------
    logging.info("Loading data from PBI-Scope...")

    sys.path.insert(0, str(Path(__file__).parent))
    from pbi import quick_connect
    from pbi.negative_examples import NegativeExampleGenerator
    from pbi_adapter import PBIAdapter
    from sklearn.model_selection import train_test_split

    retriever = quick_connect()
    adapter = PBIAdapter(
        retriever,
        bacterium_threshold=args.bacterium_threshold,
        phage_threshold=args.phage_threshold,
        bacterium_min_length=args.bacterium_min_length,
        phage_min_length=args.phage_min_length,
    )

    # Phase 1: Query pair IDs only (fast — no sequences fetched)
    logging.info("Querying pair IDs from database...")
    all_pairs = adapter.get_pair_ids_only()
    logging.info(f"Found {len(all_pairs)} pairs in database")

    # Classify pairs by interaction type (before applying limit!)
    logging.info("Classifying pairs by interaction type...")
    positive_pairs, private_negatives = adapter.classify_pairs_by_interaction(all_pairs)

    # Apply limit to positive pairs only (private negatives are always included)
    if args.limit and len(positive_pairs) > args.limit:
        logging.info(
            f"Limiting positive pairs: {len(positive_pairs)} -> {args.limit} "
            f"(private negatives always included)"
        )
        positive_pairs = positive_pairs.sample(n=args.limit, random_state=42).reset_index(drop=True)

    logging.info(f"Positive pairs: {len(positive_pairs)}")
    logging.info(f"True negatives from private data: {len(private_negatives)}")

    # Generate synthetic negatives
    logging.info(f"Generating synthetic negatives (ratio={args.negative_ratio})...")
    neg_gen = NegativeExampleGenerator(retriever)
    generated_negatives = neg_gen.generate_random_negatives(
        positive_pairs, ratio=args.negative_ratio
    )
    generated_negatives["negative_source"] = "generated"
    logging.info(f"Generated {len(generated_negatives)} synthetic negative pairs")

    # Combine all negatives
    if len(private_negatives) > 0 and len(generated_negatives) > 0:
        negative_pairs = pd.concat([private_negatives, generated_negatives], ignore_index=True)
    elif len(private_negatives) > 0:
        negative_pairs = private_negatives
    else:
        negative_pairs = generated_negatives

    logging.info(
        f"Total negatives: {len(negative_pairs)} "
        f"({len(private_negatives)} private_data + "
        f"{len(generated_negatives)} generated)"
    )

    # Prepare training data
    logging.info("Preparing training data (fetching sequences, padding, encoding)...")
    couples, labels = adapter.prepare_training_data(positive_pairs, negative_pairs)
    logging.info(
        f"Total pairs: {len(couples)} "
        f"({int(labels.sum())} positive, {int(len(labels) - labels.sum())} negative)"
    )

    # Train/validation/test split
    X_train, X_test, y_train, y_test = train_test_split(
        couples, labels, stratify=labels, test_size=0.3, shuffle=True, random_state=42
    )
    X_valid, X_test, y_valid, y_test = train_test_split(
        X_test, y_test, stratify=y_test, test_size=0.5, shuffle=True, random_state=42
    )
    logging.info(f"Split: train={len(X_train)}, valid={len(X_valid)}, test={len(X_test)}")

    # -----------------------------------------------------------------------
    # Build model
    # -----------------------------------------------------------------------
    logging.info("Building PERPHECT model...")
    model = build_model(args.bacterium_threshold, args.phage_threshold)

    import keras
    optimizer = keras.optimizers.Adam(learning_rate=args.learning_rate)
    model.compile(optimizer, "binary_crossentropy", metrics=["accuracy"], jit_compile=False)
    model.summary()

    # -----------------------------------------------------------------------
    # Training
    # -----------------------------------------------------------------------
    logging.info("Starting training...")

    # Create generators
    train_gen = adapter.create_tf_generator(
        X_train, y_train, batch_size=args.batch_size, shuffle=True
    )
    valid_gen = adapter.create_tf_generator(
        X_valid, y_valid, batch_size=args.batch_size, shuffle=False
    )
    valid_steps = math.ceil(len(X_valid) / args.batch_size)

    # Callbacks
    callbacks = [
        keras.callbacks.ModelCheckpoint(
            filepath=str(output_dir / "model_best.keras"),
            monitor="val_loss",
            save_best_only=True,
        ),
        keras.callbacks.EarlyStopping(
            monitor="val_loss",
            patience=args.patience,
            restore_best_weights=True,
        ),
        keras.callbacks.LearningRateScheduler(
            lambda epoch: step_decay(epoch, args.learning_rate)
        ),
        keras.callbacks.CSVLogger(str(output_dir / "training_log.csv")),
    ]

    # Train
    history = model.fit(
        train_gen,
        steps_per_epoch=args.steps_per_epoch,
        epochs=args.epochs,
        validation_data=valid_gen,
        validation_steps=valid_steps,
        callbacks=callbacks,
    )

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    logging.info("Saving results...")

    # Save final model
    model.save(str(output_dir / "model_final.keras"))

    # Save training plots
    try:
        from plotting_utils import historic_plots
        import matplotlib
        matplotlib.use("Agg")

        acc_fig, loss_fig = historic_plots(history)
        acc_fig.axes[0].set_title("Accuracy")
        loss_fig.axes[0].set_title("Loss")
        acc_fig.savefig(str(output_dir / "accuracy.png"), dpi=150, bbox_inches="tight")
        loss_fig.savefig(str(output_dir / "val_loss.png"), dpi=150, bbox_inches="tight")
        import matplotlib.pyplot as plt
        plt.close("all")
        logging.info("Training plots saved")
    except Exception as e:
        logging.warning(f"Could not save plots: {e}")

    # -----------------------------------------------------------------------
    # Test evaluation
    # -----------------------------------------------------------------------
    logging.info("Evaluating on test set...")

    test_gen = adapter.create_tf_generator(
        X_test, y_test, batch_size=args.batch_size, shuffle=False
    )
    test_steps = math.ceil(len(X_test) / args.batch_size)
    test_predictions = model.predict(test_gen, steps=test_steps)
    test_pred_labels = (test_predictions.flatten() > 0.5).astype(int)

    # Classification report
    from sklearn.metrics import classification_report, confusion_matrix
    report = classification_report(
        y_test, test_pred_labels, target_names=["Negative", "Positive"]
    )
    logging.info(f"\nClassification Report:\n{report}")

    # Save test results
    results_df = pd.DataFrame({
        "bacterium_id": X_test[:, 0],
        "phage_id": X_test[:, 1],
        "observations": y_test,
        "predictions": test_predictions.flatten(),
        "prediction_labels": test_pred_labels,
    })
    results_df.to_csv(str(output_dir / "results_test_set.csv"), index=False)

    # Save summary
    n_private_neg = len(private_negatives)
    n_generated_neg = len(generated_negatives)
    summary = {
        "run_name": run_name,
        "total_pairs": len(couples),
        "positive_pairs": int(labels.sum()),
        "negative_pairs": int(len(labels) - labels.sum()),
        "negative_sources": {
            "private_data": n_private_neg,
            "generated": n_generated_neg,
        },
        "train_size": len(X_train),
        "valid_size": len(X_valid),
        "test_size": len(X_test),
        "epochs_trained": len(history.history["loss"]),
        "best_val_loss": float(min(history.history["val_loss"])),
        "best_val_accuracy": float(max(history.history["val_accuracy"])),
        "test_report": report,
    }
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logging.info("=" * 60)
    logging.info("Training complete!")
    logging.info(f"Results saved to: {output_dir}")
    logging.info(f"Best model: {output_dir / 'model_best.keras'}")
    logging.info(f"Summary: {output_dir / 'summary.json'}")
    logging.info("=" * 60)


if __name__ == "__main__":
    main()
