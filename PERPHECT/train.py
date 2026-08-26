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
        logging.info(f"TensorFlow version: {tf.__version__}")
        logging.info(f"Built with CUDA: {tf.test.is_built_with_cuda()}")
        try:
            logging.info(f"cuDNN version: {tf.sysconfig.get_build_info().get('cudnn_version', 'unknown')}")
        except Exception:
            pass
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


def focal_loss(gamma=2.0, alpha=0.25):
    """Focal loss for imbalanced binary classification.

    Down-weights easy negatives automatically so the model focuses on
    hard examples.  No class_weight or sample_weight needed.

    Args:
        gamma: Focusing parameter (higher = more focus on hard examples).
               0.0 recovers standard binary crossentropy.
        alpha: Weight for the positive class.  0.25 is the standard value
               from the original paper (Lin et al., 2017).
    """
    import keras
    import tensorflow as tf

    def loss(y_true, y_pred):
        epsilon = keras.backend.epsilon()
        y_pred = tf.clip_by_value(y_pred, epsilon, 1.0 - epsilon)
        bce = -(
            y_true * tf.math.log(y_pred)
            + (1.0 - y_true) * tf.math.log(1.0 - y_pred)
        )
        p_t = y_true * y_pred + (1.0 - y_true) * (1.0 - y_pred)
        alpha_t = y_true * alpha + (1.0 - y_true) * (1.0 - alpha)
        return tf.reduce_mean(alpha_t * tf.pow(1.0 - p_t, gamma) * bce)

    return loss


def freeze_base_layers(model, up_to_layer="concatenated_features"):
    """Freeze base encoder layers up to (and including) the specified layer."""
    for layer in model.layers:
        layer.trainable = False
        if layer.name == up_to_layer:
            break
    trainable_count = sum(1 for l in model.layers if l.trainable)
    total_count = len(model.layers)
    logging.info(f"Frozen layers up to '{up_to_layer}': {total_count - trainable_count}/{total_count} frozen, {trainable_count} trainable")
    return trainable_count


def load_or_build_model(args, bacterium_threshold, phage_threshold, is_finetuning=False):
    """Load pre-trained model or build new one. Returns (model, is_finetuning_flag)."""
    import keras
    
    if args.pretrained_model:
        logging.info(f"Loading pre-trained model from {args.pretrained_model}")
        model = keras.models.load_model(
            args.pretrained_model,
            custom_objects={"focal_loss": focal_loss}
        )
        is_finetuning = True
        
        if args.freeze_base:
            freeze_base_layers(model, args.freeze_up_to)
        
        # Recompile with fine-tuning LR
        lr = args.fine_tune_lr if is_finetuning else args.learning_rate
        optimizer = keras.optimizers.Adam(learning_rate=lr)
        model.compile(
            optimizer,
            focal_loss(gamma=args.focal_gamma, alpha=args.focal_alpha),
            metrics=["accuracy", keras.metrics.AUC(name="auc")],
            jit_compile=False,
        )
        logging.info(f"Model loaded for {'fine-tuning' if is_finetuning else 'training'} with LR={lr}")
        return model, True
    else:
        model = build_model(bacterium_threshold, phage_threshold)
        optimizer = keras.optimizers.Adam(learning_rate=args.learning_rate)
        model.compile(
            optimizer,
            focal_loss(gamma=args.focal_gamma, alpha=args.focal_alpha),
            metrics=["accuracy", keras.metrics.AUC(name="auc")],
            jit_compile=False,
        )
        logging.info(f"New model built for training with LR={args.learning_rate}")
        return model, False


def _train_fold(fold_num, adapter, X_train, y_train, X_valid, y_valid,
                X_test, y_test, args, output_dir, test_pair_ids=None):
    """Train and evaluate a single fold. Returns dict of fold metrics."""
    import keras

    class ValidationProgressCallback(keras.callbacks.Callback):
        """Logs validation progress so the user knows it's not stuck."""
        def __init__(self, validation_steps):
            super().__init__()
            self.validation_steps = validation_steps
        def on_test_begin(self, logs=None):
            logging.info(f"  Validation: 0/{self.validation_steps} steps")
        def on_test_batch_end(self, batch, logs=None):
            if (batch + 1) % max(1, self.validation_steps // 10) == 0 or batch + 1 == self.validation_steps:
                logging.info(f"  Validation: {batch + 1}/{self.validation_steps} steps")

    fold_dir = output_dir / f"fold_{fold_num}"
    fold_dir.mkdir(parents=True, exist_ok=True)

    # Build or load model
    model, is_finetuning_fold = load_or_build_model(
        args, args.bacterium_threshold, args.phage_threshold, is_finetuning=False
    )

    # Generators
    train_gen = adapter.create_tf_generator(
        X_train, y_train, batch_size=args.batch_size, shuffle=True
    )
    valid_gen = adapter.create_tf_generator(
        X_valid, y_valid, batch_size=args.batch_size, shuffle=False
    )
    valid_steps = math.ceil(len(X_valid) / args.batch_size)

    # Auto-calculate steps_per_epoch if not set
    steps_per_epoch = args.steps_per_epoch
    if steps_per_epoch is None:
        steps_per_epoch = math.ceil(len(X_train) / args.batch_size)

    # Callbacks
    callbacks = [
        keras.callbacks.ModelCheckpoint(
            filepath=str(fold_dir / "model_best.keras"),
            monitor="val_auc", mode="max", save_best_only=True,
        ),
        keras.callbacks.EarlyStopping(
            monitor="val_auc", mode="max", patience=args.patience, restore_best_weights=True,
        ),
        keras.callbacks.LearningRateScheduler(
            lambda epoch: step_decay(epoch, args.learning_rate)
        ),
        keras.callbacks.CSVLogger(str(fold_dir / "training_log.csv")),
        ValidationProgressCallback(valid_steps),
    ]

    # Train
    history = model.fit(
        train_gen, steps_per_epoch=steps_per_epoch,
        epochs=args.epochs, validation_data=valid_gen,
        validation_steps=valid_steps, callbacks=callbacks,
    )

    # Save final model
    model.save(str(fold_dir / "model_final.keras"))

    # Test evaluation
    test_gen = adapter.create_tf_generator(
        X_test, y_test, batch_size=args.batch_size, shuffle=False
    )
    test_steps = math.ceil(len(X_test) / args.batch_size)
    test_predictions = model.predict(test_gen, steps=test_steps)
    test_pred_labels = (test_predictions.flatten() > 0.5).astype(int)

    from sklearn.metrics import classification_report, matthews_corrcoef, confusion_matrix
    report = classification_report(
        y_test, test_pred_labels, target_names=["Negative", "Positive"]
    )
    mcc = matthews_corrcoef(y_test, test_pred_labels)
    tn, fp, fn, tp = confusion_matrix(y_test, test_pred_labels).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    logging.info(f"  MCC: {mcc:.4f}  Sensitivity: {sensitivity:.4f}  Specificity: {specificity:.4f}")

    results_data = {
        "bacterium_id": X_test[:, 0],
        "phage_id_int": X_test[:, 1],
        "observations": y_test,
        "predictions": test_predictions.flatten(),
        "prediction_labels": test_pred_labels,
    }
    if test_pair_ids is not None:
        results_data["host_id"] = test_pair_ids["host_id"].values
        results_data["phage_id"] = test_pair_ids["phage_id"].values
    results_df = pd.DataFrame(results_data)
    results_df.to_csv(str(fold_dir / "results_test_set.csv"), index=False)

    # Save summary
    summary = {
        "fold": fold_num,
        "train_size": len(X_train),
        "valid_size": len(X_valid),
        "test_size": len(X_test),
        "epochs_trained": len(history.history["loss"]),
        "best_val_loss": float(min(history.history["val_loss"])),
        "best_val_auc": float(max(history.history["val_auc"])),
        "best_val_accuracy": float(max(history.history["val_accuracy"])),
        "test_mcc": float(mcc),
        "test_sensitivity": float(sensitivity),
        "test_specificity": float(specificity),
        "test_report": report,
    }
    with open(fold_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logging.info(f"Fold {fold_num}: val_auc={summary['best_val_auc']:.4f}, "
                 f"val_loss={summary['best_val_loss']:.4f}, "
                 f"test_mcc={mcc:.4f}, "
                 f"epochs={summary['epochs_trained']}")

    return summary


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
    parser.add_argument("--steps-per-epoch", type=int, default=None, help="Batches per epoch (None = cover full training set)")
    parser.add_argument("--patience", type=int, default=5, help="Early stopping patience (default: 5)")
    parser.add_argument("--learning-rate", type=float, default=0.0004, help="Initial learning rate (default: 0.0004)")

    # Data parameters
    parser.add_argument("--limit", type=int, default=None, help="Limit positive pairs (None = all)")
    parser.add_argument("--negative-ratio", type=float, default=1.0,
                        help="Max synthetic negatives as multiple of private negatives "
                             "(1.0 = match private count, 2.0 = 2x, 0.0 = no synthetic). "
                             "Default: 1.0 (gives ~50/50 balance before filtering)")
    parser.add_argument("--focal-alpha", type=float, default=0.25,
                        help="Focal loss positive-class weight (default: 0.25)")
    parser.add_argument("--focal-gamma", type=float, default=2.0,
                        help="Focal loss focusing parameter (default: 2.0)")
    parser.add_argument("--bacterium-threshold", type=int, default=7_000_000, help="Max bacteria sequence length")
    parser.add_argument("--phage-threshold", type=int, default=200_000, help="Max phage sequence length")
    parser.add_argument("--bacterium-min-length", type=int, default=150_000, help="Min bacteria length to keep")
    parser.add_argument("--phage-min-length", type=int, default=1_500, help="Min phage length to keep")

    # Output parameters
    parser.add_argument("--output-dir", type=str, default="/results", help="Output directory (default: /results, mapped to ./outputs on host)")
    parser.add_argument("--run-name", type=str, default=None, help="Run name (default: timestamp)")

    # Configuration
    parser.add_argument("--config", type=str, default=None, help="YAML config file")
    parser.add_argument("--profile", type=str, default=None, help="Profile name from config (e.g., 'pretrain', 'finetune')")
    parser.add_argument("--cross-validate", type=int, default=None, help="K folds for stratified K-fold CV (overrides config)")
    
    # Pre-trained model / fine-tuning
    parser.add_argument("--pretrained-model", type=str, default=None, help="Path to pre-trained .keras model to load (enables fine-tuning mode)")
    parser.add_argument("--freeze-base", action="store_true", help="Freeze CNN encoder layers, only train classification head")
    parser.add_argument("--freeze-up-to", type=str, default="concatenated_features", help="Freeze layers up to this layer name (default: concatenated_features)")
    parser.add_argument("--fine-tune-lr", type=float, default=0.0001, help="Learning rate for fine-tuning (default: 0.0001)")
    parser.add_argument("--fine-tune-epochs", type=int, default=5, help="Epochs for fine-tuning (default: 5)")
    parser.add_argument("--finetuned-model-name", type=str, default="model_finetuned_best.keras", help="Output name for fine-tuned best model")
    
    parser.add_argument("--no-gpu", action="store_true", help="Force CPU even if GPU available")
    parser.add_argument("--gpu-device", type=int, default=0, help="GPU device index (default: 0, use 1 for second GPU)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument("--log-file", type=str, default=None, help="Log to file in addition to console")

    args = parser.parse_args()

    # Load config file if provided
    if args.config:
        import yaml
        with open(args.config) as f:
            config = yaml.safe_load(f)

        # Resolve which config section to use: defaults + optional profile overlay
        if "defaults" in config:
            merged = config["defaults"]
        else:
            merged = config  # backward compat: flat config treated as defaults

        if args.profile:
            profiles = config.get("profiles", {})
            if args.profile not in profiles:
                available = list(profiles.keys()) if profiles else "(none)"
                raise ValueError(
                    f"Profile '{args.profile}' not found in {args.config}. "
                    f"Available profiles: {available}"
                )
            profile = profiles[args.profile]
            for section in ["training", "data", "output", "gpu"]:
                if section in profile:
                    if section not in merged:
                        merged[section] = {}
                    merged[section].update(profile[section])
            logging.info(f"Using profile: {args.profile}")

        # Apply config values to args (only when CLI arg was not explicitly set)
        for section in ["training", "data", "output", "gpu"]:
            if section in merged:
                for key, value in merged[section].items():
                    attr = key.replace("-", "_")
                    current = getattr(args, attr, None)
                    if current is None:
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
    from sklearn.model_selection import train_test_split, StratifiedKFold

    retriever = quick_connect()
    adapter = PBIAdapter(
        retriever,
        bacterium_threshold=args.bacterium_threshold,
        phage_threshold=args.phage_threshold,
        bacterium_min_length=args.bacterium_min_length,
        phage_min_length=args.phage_min_length,
    )

    # Parse exclude-sources (list from config or comma-separated string)
    exclude_sources = None
    if args.exclude_sources:
        if isinstance(args.exclude_sources, list):
            exclude_sources = [s.strip() for s in args.exclude_sources if s.strip()]
        else:
            exclude_sources = [s.strip() for s in args.exclude_sources.split(",") if s.strip()]
        logging.info(f"Excluding sources: {exclude_sources}")

    is_finetuning = args.pretrained_model is not None

    # Phase 1: Query pair IDs only (fast — no sequences fetched)
    logging.info("Querying pair IDs from database...")
    
    if is_finetuning:
        # FINE-TUNING MODE: Only use the excluded sources
        if not exclude_sources:
            raise ValueError("--exclude-sources is required for fine-tuning mode")
        
        # Query ONLY the specified sources for fine-tuning
        placeholders = ", ".join(["?" for _ in exclude_sources])
        query = f"""
        SELECT DISTINCT pha.Phage_ID, pha.Host_ID
        FROM phage_host_associations pha
        JOIN fact_phages p ON pha.Phage_ID = p.Phage_ID
        WHERE p.Source_DB IN ({placeholders})
        """
        query += " ORDER BY MD5(Phage_ID || Host_ID)"
        all_pairs = retriever.conn.execute(query, exclude_sources).fetchdf()
        logging.info(f"Fine-tuning mode: Found {len(all_pairs)} pairs from sources {exclude_sources}")
    else:
        # PRE-TRAINING MODE: Exclude specified sources
        all_pairs = adapter.get_pair_ids_only(shuffle=True, exclude_sources=exclude_sources)
        logging.info(f"Pre-training mode: Found {len(all_pairs)} pairs in database (excluded sources: {exclude_sources})")

    # Exclude held-out test pairs if specified
    if args.exclude_ids:
        exclude_path = Path(args.exclude_ids)
        if exclude_path.exists():
            exclude_df = pd.read_csv(exclude_path)
            exclude_set = set(zip(exclude_df["Phage_ID"], exclude_df["Host_ID"]))
            before = len(all_pairs)
            all_pairs = all_pairs[
                ~all_pairs.apply(
                    lambda r: (r["Phage_ID"], r["Host_ID"]) in exclude_set, axis=1
                )
            ].reset_index(drop=True)
            logging.info(
                f"Excluded {before - len(all_pairs)} test pairs from {exclude_path} "
                f"({len(all_pairs)} remaining)"
            )
        else:
            raise FileNotFoundError(
                f"--exclude-ids file not found: {exclude_path}. "
                "Run 01_prepare_test_set.ipynb first to create the excluded pairs CSV."
            )

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
    logging.info(
        f"Existing negatives from private data: {len(private_negatives)} "
        f"(pairs explicitly marked as 'no interaction' in private_interactions)"
    )

    # Generate synthetic negatives (capped by --negative-ratio × private negatives)
    n_private = len(private_negatives)
    if n_private > 0:
        max_generated = int(n_private * args.negative_ratio)
        logging.info(
            f"Generating synthetic negatives "
            f"(capped at {args.negative_ratio:.1f}x existing negatives = {max_generated})"
        )
    else:
        max_generated = len(positive_pairs)
        logging.info(
            f"No existing negatives in private data — generating synthetic negatives "
            f"to create a negative class"
        )
    neg_gen = NegativeExampleGenerator(retriever)
    target_count = min(len(positive_pairs), max_generated)
    effective_ratio = target_count / len(positive_pairs) if len(positive_pairs) > 0 else 0
    logging.info(f"  Target: {target_count} synthetic negatives (ratio={effective_ratio:.3f})")
    generated_negatives = neg_gen.generate_random_negatives(
        positive_pairs, ratio=effective_ratio
    )

    # Deduplicate against private negatives
    # NOTE: Generated negatives are NOT deduped against the --exclude-ids test set.
    # Collision probability is negligible (billions of possible phage-host pairs).
    if len(private_negatives) > 0:
        private_neg_set = set(
            zip(private_negatives["Phage_ID"], private_negatives["Host_ID"])
        )
        before_dedup = len(generated_negatives)
        generated_negatives = generated_negatives[
            ~generated_negatives.apply(
                lambda r: (r["Phage_ID"], r["Host_ID"]) in private_neg_set, axis=1
            )
        ].reset_index(drop=True)
        if before_dedup - len(generated_negatives) > 0:
            logging.info(
                f"  Removed {before_dedup - len(generated_negatives)} duplicates "
                f"against existing private negatives"
            )

    generated_negatives["negative_source"] = "generated"
    logging.info(f"  Generated: {len(generated_negatives)} synthetic negative pairs")

    # Combine all negatives
    if len(private_negatives) > 0 and len(generated_negatives) > 0:
        negative_pairs = pd.concat([private_negatives, generated_negatives], ignore_index=True)
    elif len(private_negatives) > 0:
        negative_pairs = private_negatives
    else:
        negative_pairs = generated_negatives

    logging.info(
        f"Final negative pool: {len(negative_pairs)} total "
        f"({len(private_negatives)} from private data + "
        f"{len(generated_negatives)} generated)"
    )

    # Prepare training data
    logging.info("Preparing training data (fetching sequences, padding, encoding)...")
    couples, labels, sources = adapter.prepare_training_data(positive_pairs, negative_pairs)
    logging.info(
        f"Total pairs: {len(couples)} "
        f"({int(labels.sum())} positive, {int(len(labels) - labels.sum())} negative"
        f")"
    )

    # Build pair_ids DataFrame for traceability (maps integer indices back to string IDs)
    reverse_host = {v: k for k, v in adapter.host_id_map.items()}
    reverse_phage = {v: k for k, v in adapter.phage_id_map.items()}
    pair_ids = pd.DataFrame({
        "host_id": [reverse_host.get(int(c[0]), str(c[0])) for c in couples],
        "phage_id": [reverse_phage.get(int(c[1]), str(c[1])) for c in couples],
        "label": labels.astype(int),
        "source": sources,
    })
    pair_ids.to_csv(str(output_dir / "pairs_all.csv"), index=False)
    logging.info(f"Saved all pair IDs to {output_dir / 'pairs_all.csv'}")

    # -----------------------------------------------------------------------
    # Train/validation/test split (stratified by label AND source)
    # -----------------------------------------------------------------------
    logging.info("Splitting data into train/validation/test (stratified)...")
    stratify_key = np.array([
        "pos" if l == 1 else f"neg_{s}"
        for l, s in zip(labels, sources)
    ])

    # -----------------------------------------------------------------------
    # Branch: Cross-validation or standard split
    # -----------------------------------------------------------------------
    if (args.cross_validate or 0) > 0:
        if args.pretrained_model:
            raise ValueError("Cross-validation (--cross-validate) is not supported with fine-tuning (--pretrained-model). Use standard split for fine-tuning.")
        # ---- K-fold cross-validation ----
        k = args.cross_validate or 0
        logging.info(f"Starting {k}-fold stratified cross-validation...")

        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
        fold_summaries = []

        for fold_num, (train_idx, test_idx) in enumerate(skf.split(couples, stratify_key), 1):
            logging.info(f"\n{'='*60}")
            logging.info(f"Fold {fold_num}/{k}")
            logging.info(f"{'='*60}")

            X_fold_train, X_fold_test = couples[train_idx], couples[test_idx]
            y_fold_train, y_fold_test = labels[train_idx], labels[test_idx]

            # Split train into train/val (80/20 of fold train)
            fold_strat = np.array([
                "pos" if l == 1 else f"neg_{s}"
                for l, s in zip(y_fold_train, sources[train_idx])
            ])
            X_fold_train, X_fold_val, y_fold_train, y_fold_val = train_test_split(
                X_fold_train, y_fold_train,
                stratify=fold_strat, test_size=0.2, shuffle=True, random_state=42
            )

            # Propagate pair_ids through this fold's splits
            fold_pair_ids = pair_ids.iloc[test_idx].reset_index(drop=True)
            fold_train_pair_ids, fold_val_pair_ids = train_test_split(
                pair_ids.iloc[train_idx].reset_index(drop=True),
                stratify=fold_strat, test_size=0.2, shuffle=True, random_state=42
            )

            summary = _train_fold(
                fold_num, adapter,
                X_fold_train, y_fold_train,
                X_fold_val, y_fold_val,
                X_fold_test, y_fold_test,
                args, output_dir,
                test_pair_ids=fold_pair_ids,
            )

            # Save per-fold pair IDs for traceability
            fold_dir = output_dir / f"fold_{fold_num}"
            fold_train_pair_ids.to_csv(str(fold_dir / "pairs_train.csv"), index=False)
            fold_val_pair_ids.to_csv(str(fold_dir / "pairs_val.csv"), index=False)
            fold_pair_ids.to_csv(str(fold_dir / "pairs_test.csv"), index=False)

            fold_summaries.append(summary)

        # Aggregate CV results
        cv_summary = {
            "run_name": run_name,
            "n_folds": k,
            "total_pairs": len(couples),
            "positive_pairs": int(labels.sum()),
            "negative_pairs": int(len(labels) - labels.sum()),
            "negative_sources": {
                "private_data": len(private_negatives),
                "generated": len(generated_negatives),
            },
            "pair_ids_file": "pairs_all.csv",
            "mean_val_auc": float(np.mean([s["best_val_auc"] for s in fold_summaries])),
            "std_val_auc": float(np.std([s["best_val_auc"] for s in fold_summaries])),
            "mean_val_loss": float(np.mean([s["best_val_loss"] for s in fold_summaries])),
            "std_val_loss": float(np.std([s["best_val_loss"] for s in fold_summaries])),
            "mean_test_mcc": float(np.mean([s["test_mcc"] for s in fold_summaries])),
            "std_test_mcc": float(np.std([s["test_mcc"] for s in fold_summaries])),
            "mean_epochs": float(np.mean([s["epochs_trained"] for s in fold_summaries])),
            "folds": fold_summaries,
        }
        with open(output_dir / "cv_summary.json", "w") as f:
            json.dump(cv_summary, f, indent=2)

        logging.info(f"\n{'='*60}")
        logging.info("Cross-validation complete!")
        logging.info(f"Val AUC:   {cv_summary['mean_val_auc']:.4f} +/- {cv_summary['std_val_auc']:.4f}")
        logging.info(f"Val loss:  {cv_summary['mean_val_loss']:.4f} +/- {cv_summary['std_val_loss']:.4f}")
        logging.info(f"Test MCC:  {cv_summary['mean_test_mcc']:.4f} +/- {cv_summary['std_test_mcc']:.4f}")
        logging.info(f"Mean epochs: {cv_summary['mean_epochs']:.1f}")
        logging.info(f"Results saved to: {output_dir}")
        logging.info(f"Summary: {output_dir / 'cv_summary.json'}")
        logging.info("=" * 60)

    else:
        # ---- Standard train/val/test split ----
        X_train, X_test, y_train, y_test, s_train, s_test = train_test_split(
            couples, labels, sources,
            stratify=stratify_key, test_size=0.3, shuffle=True, random_state=42
        )
        pairs_train, pairs_test = train_test_split(
            pair_ids, stratify=stratify_key, test_size=0.3, shuffle=True, random_state=42
        )

        stratify_key_test = np.array([
            "pos" if l == 1 else f"neg_{s}"
            for l, s in zip(y_test, s_test)
        ])
        X_valid, X_test, y_valid, y_test, s_valid, s_test = train_test_split(
            X_test, y_test, s_test,
            stratify=stratify_key_test, test_size=0.5, shuffle=True, random_state=42
        )
        pairs_val, pairs_test = train_test_split(
            pairs_test, stratify=stratify_key_test, test_size=0.5, shuffle=True, random_state=42
        )

        # Save per-split pair IDs for traceability
        pairs_train.to_csv(str(output_dir / "pairs_train.csv"), index=False)
        pairs_val.to_csv(str(output_dir / "pairs_val.csv"), index=False)
        pairs_test.to_csv(str(output_dir / "pairs_test.csv"), index=False)

        logging.info(f"Split: train={len(X_train)}, valid={len(X_valid)}, test={len(X_test)}")

        for split_name, split_sources in [("train", s_train), ("valid", s_valid), ("test", s_test)]:
            unique, counts = np.unique(split_sources, return_counts=True)
            dist = dict(zip(unique, counts))
            logging.info(f"  {split_name}: {dist}")

        # Build model
        logging.info("Building/loading PERPHECT model...")
        model, is_finetuning = load_or_build_model(
            args, args.bacterium_threshold, args.phage_threshold, is_finetuning=args.pretrained_model is not None
        )
        model.summary()

        class ValidationProgressCallback(keras.callbacks.Callback):
            """Logs validation progress so the user knows it's not stuck."""
            def __init__(self, validation_steps):
                super().__init__()
                self.validation_steps = validation_steps
            def on_test_begin(self, logs=None):
                logging.info(f"  Validation: 0/{self.validation_steps} steps")
            def on_test_batch_end(self, batch, logs=None):
                if (batch + 1) % max(1, self.validation_steps // 10) == 0 or batch + 1 == self.validation_steps:
                    logging.info(f"  Validation: {batch + 1}/{self.validation_steps} steps")

        # Training
        logging.info("Starting training...")

        train_gen = adapter.create_tf_generator(
            X_train, y_train, batch_size=args.batch_size, shuffle=True
        )
        valid_gen = adapter.create_tf_generator(
            X_valid, y_valid, batch_size=args.batch_size, shuffle=False
        )
        valid_steps = math.ceil(len(X_valid) / args.batch_size)

        # Determine training parameters based on mode
        if is_finetuning:
            train_epochs = args.fine_tune_epochs
            train_lr = args.fine_tune_lr
            best_model_name = args.finetuned_model_name
            logging.info(f"Fine-tuning mode: epochs={train_epochs}, lr={train_lr}, best_model={best_model_name}")
        else:
            train_epochs = args.epochs
            train_lr = args.learning_rate
            best_model_name = "model_best.keras"
            logging.info(f"Training mode: epochs={train_epochs}, lr={train_lr}")

        # Auto-calculate steps_per_epoch if not set
        steps_per_epoch = args.steps_per_epoch
        if steps_per_epoch is None:
            steps_per_epoch = math.ceil(len(X_train) / args.batch_size)
            logging.info(f"Auto steps_per_epoch: {steps_per_epoch} (full training set)")

        callbacks = [
            keras.callbacks.ModelCheckpoint(
                filepath=str(output_dir / best_model_name),
                monitor="val_auc", mode="max", save_best_only=True,
            ),
            keras.callbacks.EarlyStopping(
                monitor="val_auc", mode="max", patience=args.patience, restore_best_weights=True,
            ),
            keras.callbacks.LearningRateScheduler(
                lambda epoch: step_decay(epoch, train_lr)
            ),
            keras.callbacks.CSVLogger(str(output_dir / "training_log.csv")),
            ValidationProgressCallback(valid_steps),
        ]

        history = model.fit(
            train_gen, steps_per_epoch=steps_per_epoch,
            epochs=train_epochs, validation_data=valid_gen,
            validation_steps=valid_steps, callbacks=callbacks,
        )

        # Save results
        logging.info("Saving results...")
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

        # Test evaluation
        logging.info("Evaluating on test set...")
        test_gen = adapter.create_tf_generator(
            X_test, y_test, batch_size=args.batch_size, shuffle=False
        )
        test_steps = math.ceil(len(X_test) / args.batch_size)
        test_predictions = model.predict(test_gen, steps=test_steps)
        test_pred_labels = (test_predictions.flatten() > 0.5).astype(int)

        from sklearn.metrics import classification_report, matthews_corrcoef, confusion_matrix
        report = classification_report(
            y_test, test_pred_labels, target_names=["Negative", "Positive"]
        )
        mcc = matthews_corrcoef(y_test, test_pred_labels)
        tn, fp, fn, tp = confusion_matrix(y_test, test_pred_labels).ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        logging.info(f"\nClassification Report:\n{report}")
        logging.info(f"MCC: {mcc:.4f}  Sensitivity: {sensitivity:.4f}  Specificity: {specificity:.4f}")

        results_df = pd.DataFrame({
            "host_id": pairs_test["host_id"].values,
            "phage_id": pairs_test["phage_id"].values,
            "bacterium_id": X_test[:, 0],
            "phage_id_int": X_test[:, 1],
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
            "pair_ids_file": "pairs_all.csv",
            "train_size": len(X_train),
            "valid_size": len(X_valid),
            "test_size": len(X_test),
            "split_sources": {
                split: {str(k): int(v) for k, v in zip(*np.unique(arr, return_counts=True))}
                for split, arr in [("train", s_train), ("valid", s_valid), ("test", s_test)]
            },
            "epochs_trained": len(history.history["loss"]),
            "best_val_loss": float(min(history.history["val_loss"])),
            "best_val_auc": float(max(history.history["val_auc"])),
            "best_val_accuracy": float(max(history.history["val_accuracy"])),
            "test_mcc": float(mcc),
            "test_sensitivity": float(sensitivity),
            "test_specificity": float(specificity),
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
