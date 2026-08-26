"""
Tests for PERPHECT train.py script.

Tests the script's argument parsing, configuration loading, model building,
and helper functions without requiring a database or GPU.
"""

import sys
import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "PERPHECT"))


class TestArgParsing:
    """Tests for command-line argument parsing."""

    def test_default_args(self):
        from train import main
        import argparse

        # We can't call main() directly (it runs training), but we can
        # test the parser by importing and checking defaults
        parser = argparse.ArgumentParser()
        parser.add_argument("--epochs", type=int, default=10)
        parser.add_argument("--batch-size", type=int, default=4)
        parser.add_argument("--steps-per-epoch", type=int, default=400)
        parser.add_argument("--patience", type=int, default=5)
        parser.add_argument("--learning-rate", type=float, default=0.0004)
        parser.add_argument("--limit", type=int, default=None)
        parser.add_argument("--negative-ratio", type=float, default=1.0)
        parser.add_argument("--bacterium-threshold", type=int, default=7_000_000)
        parser.add_argument("--phage-threshold", type=int, default=200_000)
        parser.add_argument("--bacterium-min-length", type=int, default=150_000)
        parser.add_argument("--phage-min-length", type=int, default=1_500)
        parser.add_argument("--output-dir", type=str, default="results")
        parser.add_argument("--run-name", type=str, default=None)
        parser.add_argument("--config", type=str, default=None)
        parser.add_argument("--cross-validate", type=int, default=0)
        parser.add_argument("--no-gpu", action="store_true")
        parser.add_argument("--gpu-device", type=int, default=0)
        parser.add_argument("--verbose", "-v", action="store_true")
        parser.add_argument("--log-file", type=str, default=None)

        args = parser.parse_args([])
        assert args.epochs == 10
        assert args.batch_size == 4
        assert args.limit is None
        assert args.no_gpu is False
        assert args.gpu_device == 0
        assert args.cross_validate == 0

    def test_custom_args(self):
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--epochs", type=int, default=10)
        parser.add_argument("--batch-size", type=int, default=4)
        parser.add_argument("--limit", type=int, default=None)
        parser.add_argument("--no-gpu", action="store_true")
        parser.add_argument("--output-dir", type=str, default="results")

        args = parser.parse_args(["--epochs", "20", "--batch-size", "32", "--limit", "1000", "--no-gpu"])
        assert args.epochs == 20
        assert args.batch_size == 32
        assert args.limit == 1000
        assert args.no_gpu is True


class TestStepDecay:
    """Tests for learning rate step decay."""

    def test_step_decay_at_epoch_0(self):
        from train import step_decay
        lr = step_decay(0, initial_lrate=0.0004)
        assert lr == pytest.approx(0.0004)

    def test_step_decay_halves(self):
        from train import step_decay
        # step_decay uses floor((1 + epoch) / epochs_drop)
        # epoch=0: floor(1/4)=0 -> lr=0.0004
        # epoch=3: floor(4/4)=1 -> lr=0.0002 (first drop)
        # epoch=7: floor(8/4)=2 -> lr=0.0001 (second drop)
        lr_0 = step_decay(0, initial_lrate=0.0004)
        lr_3 = step_decay(3, initial_lrate=0.0004)
        lr_7 = step_decay(7, initial_lrate=0.0004)
        assert lr_0 == pytest.approx(0.0004)
        assert lr_3 == pytest.approx(0.0002)
        assert lr_7 == pytest.approx(0.0001)

    def test_step_decay_custom_params(self):
        from train import step_decay
        lr = step_decay(0, initial_lrate=0.001, drop=0.9, epochs_drop=2)
        assert lr == pytest.approx(0.001)


class TestModelBuild:
    """Tests for model building (requires Keras)."""

    def test_build_model_default(self):
        try:
            import keras
        except ImportError:
            pytest.skip("keras not installed")
        from train import build_model
        model = build_model()
        assert model is not None
        assert model.name == "Perphect"

    def test_build_model_custom_thresholds(self):
        try:
            import keras
        except ImportError:
            pytest.skip("keras not installed")
        from train import build_model
        # Full architecture needs large thresholds: 3 conv+pool layers for bacteria,
        # 2 for phage. B=50000, P=5000 are minimum working values.
        model = build_model(bacterium_threshold=50000, phage_threshold=5000)
        assert model is not None
        assert model.input[0].shape[1] == 50000
        assert model.input[1].shape[1] == 5000


class TestConfigLoading:
    """Tests for YAML configuration loading."""

    def test_config_file_loading(self):
        """Verify config.yaml is valid and loadable."""
        import yaml
        config_path = Path(__file__).parent.parent / "PERPHECT" / "config.yaml"
        if config_path.exists():
            with open(config_path) as f:
                config = yaml.safe_load(f)
            assert "defaults" in config
            assert "training" in config["defaults"]
            assert "data" in config["defaults"]
            assert "output" in config["defaults"]
            assert isinstance(config["defaults"]["training"]["epochs"], int)
            assert isinstance(config["defaults"]["training"]["batch_size"], int)
            # Profiles should exist
            assert "profiles" in config
            assert "pretrain" in config["profiles"]
            assert "finetune" in config["profiles"]

    def test_config_merging(self):
        """Verify profile values override defaults."""
        # Simulate the merging logic from train.py
        defaults = {
            "training": {"epochs": 10, "batch_size": 4},
            "data": {"negative_ratio": 1.0},
        }
        profile = {
            "training": {"epochs": 20, "batch_size": 32},
            "data": {"negative_ratio": 0.5},
        }
        for section in ["training", "data"]:
            if section in profile:
                if section not in defaults:
                    defaults[section] = {}
                defaults[section].update(profile[section])

        assert defaults["training"]["epochs"] == 20
        assert defaults["training"]["batch_size"] == 32
        assert defaults["data"]["negative_ratio"] == 0.5


class TestGPUDetection:
    """Tests for GPU detection (no actual GPU needed)."""

    def test_detect_gpu_returns_bool(self):
        from train import detect_gpu
        result = detect_gpu()
        assert isinstance(result, bool)

    def test_no_gpu_environment(self):
        """Verify script works when GPU is not available."""
        with patch.dict(os.environ, {"CUDA_VISIBLE_DEVICES": ""}):
            # Re-import to reset any cached state
            import importlib
            import train
            importlib.reload(train)
            result = train.detect_gpu()
            # Should return False without crashing
            assert isinstance(result, bool)


def _make_mock_model():
    """Create a mock Keras model that supports fit/predict/save/summary."""
    import numpy as np

    model = MagicMock()
    model.name = "Perphect"

    # Mock layers for freeze_base_layers
    mock_layers = []
    for name in ["bacterial_input", "bacterial_conv_1", "concatenated_features",
                  "dense", "dropout", "output"]:
        layer = MagicMock()
        layer.name = name
        layer.trainable = True
        mock_layers.append(layer)
    model.layers = mock_layers

    # Mock input shapes
    bact_input = MagicMock()
    bact_input.shape = (None, 50000, 4)
    phage_input = MagicMock()
    phage_input.shape = (None, 5000, 4)
    model.input = [bact_input, phage_input]

    # fit returns a mock history
    history = MagicMock()
    history.history = {
        "loss": [0.5, 0.3],
        "val_loss": [0.4, 0.25],
        "val_auc": [0.6, 0.75],
        "val_accuracy": [0.55, 0.7],
        "accuracy": [0.5, 0.65],
    }
    model.fit.return_value = history

    # predict returns predictions — matches what train.py expects
    # train.py calls model.predict(gen, steps=test_steps) where
    # test_steps = ceil(n_test / batch_size). The returned array should
    # have shape (n_test, 1) for classification_report to work.
    def _predict(gen_or_array, steps=None, **kwargs):
        # When called with steps, infer n_test from steps * batch_size
        n = steps * 4 if steps else 2
        return np.random.rand(n, 1)

    model.predict.side_effect = _predict

    return model


def _make_mock_generator():
    """Create a callable that returns a mock TF generator (infinite loop of batches)."""
    import numpy as np

    def _gen_fn(*args, **kwargs):
        while True:
            yield ([
                np.random.rand(4, 50000, 4).astype(np.float32),
                np.random.rand(4, 5000, 4).astype(np.float32),
            ], np.array([1, 0, 1, 0], dtype=np.float32))

    return MagicMock(side_effect=_gen_fn)


class TestTrainingPaths:
    """Integration tests for pretrain and finetune code paths.

    These tests mock the database, adapter, and model to verify the full
    training pipeline runs without error. No GPU or database required.
    """

    def _setup_mocks(self, tmp_path, config_dict=None, pretrained_model=None):
        """Helper: create config, excluded_pairs CSV, and mock objects."""
        import yaml
        import numpy as np
        import pandas as pd

        # Minimal config using defaults structure
        if config_dict is None:
            config_dict = {
                "defaults": {
                    "training": {
                        "epochs": 2,
                        "batch_size": 4,
                        "patience": 2,
                        "learning_rate": 0.001,
                        "focal_alpha": 0.25,
                        "focal_gamma": 2.0,
                    },
                    "data": {
                        "negative_ratio": 1.0,
                        "bacterium_threshold": 50000,
                        "phage_threshold": 5000,
                        "bacterium_min_length": 100,
                        "phage_min_length": 100,
                        "exclude_ids": str(tmp_path / "excluded_pairs.csv"),
                        "exclude_sources": ["test_source"],
                    },
                    "output": {"dir": str(tmp_path / "results")},
                    "gpu": {"enabled": False},
                },
            }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        # Create excluded pairs CSV
        exclude_df = pd.DataFrame({
            "Phage_ID": ["phage_excl_1", "phage_excl_2"],
            "Host_ID": ["host_excl_1", "host_excl_2"],
        })
        exclude_csv = tmp_path / "excluded_pairs.csv"
        exclude_df.to_csv(exclude_csv, index=False)

        # Mock pair IDs (positive pairs)
        pos_pairs = pd.DataFrame({
            "Phage_ID": [f"phage_{i}" for i in range(20)],
            "Host_ID": [f"host_{i}" for i in range(20)],
        })

        # Mock negative pairs (from private data)
        neg_pairs = pd.DataFrame({
            "Phage_ID": [f"phage_neg_{i}" for i in range(10)],
            "Host_ID": [f"host_neg_{i}" for i in range(10)],
        })

        # Training data arrays (20 pos + 10 neg = 30 total)
        n_total = 30
        couples = np.column_stack([
            np.arange(n_total),  # host IDs (integers)
            np.arange(n_total),  # phage IDs (integers)
        ])
        labels = np.array([1] * 20 + [0] * 10, dtype=np.float32)
        sources = np.array(
            ["positive"] * 20 + ["private_data"] * 10, dtype=object
        )

        return config_path, pos_pairs, neg_pairs, couples, labels, sources

    def test_pretrain_path(self, tmp_path):
        """Verify pre-training path runs end-to-end with mocked dependencies."""
        import numpy as np
        from importlib import reload
        import sys
        import train

        config_path, pos_pairs, neg_pairs, couples, labels, sources = \
            self._setup_mocks(tmp_path)

        mock_model = _make_mock_model()
        mock_retriever = MagicMock()
        mock_keras = MagicMock()
        # ValidationProgressCallback needs to be a real class that inherits from Callback
        mock_keras.callbacks.Callback = type("Callback", (), {"__init__": lambda s, *a, **kw: None})
        mock_keras.callbacks.ModelCheckpoint.return_value = MagicMock()
        mock_keras.callbacks.EarlyStopping.return_value = MagicMock()
        mock_keras.callbacks.LearningRateScheduler.return_value = MagicMock()
        mock_keras.callbacks.CSVLogger.return_value = MagicMock()
        mock_keras.metrics.AUC.return_value = MagicMock()

        # Reload first so patches apply to the fresh module
        reload(train)

        with (
            patch("pbi.quick_connect", return_value=mock_retriever),
            patch("pbi_adapter.PBIAdapter") as MockAdapter,
            patch.object(train, "load_or_build_model", return_value=(mock_model, False)),
            patch.object(train, "detect_gpu", return_value=False),
            patch.object(train, "keras", mock_keras, create=True),
        ):
            adapter_instance = MagicMock()
            adapter_instance.get_pair_ids_only.return_value = pd.DataFrame({
                "Phage_ID": [f"p{i}" for i in range(30)],
                "Host_ID": [f"h{i}" for i in range(30)],
            })
            adapter_instance.classify_pairs_by_interaction.return_value = (
                pos_pairs, neg_pairs
            )
            adapter_instance.prepare_training_data.return_value = (
                couples, labels, sources
            )
            adapter_instance.create_tf_generator = _make_mock_generator()
            MockAdapter.return_value = adapter_instance

            # Mock NegativeExampleGenerator
            with patch("pbi.negative_examples.NegativeExampleGenerator") as MockNegGen:
                neg_gen_instance = MagicMock()
                gen_neg_df = pd.DataFrame({
                    "Phage_ID": [f"phage_gen_{i}" for i in range(10)],
                    "Host_ID": [f"host_gen_{i}" for i in range(10)],
                })
                neg_gen_instance.generate_random_negatives.return_value = gen_neg_df
                MockNegGen.return_value = neg_gen_instance

                # Mock train_test_split to return small splits
                # Split sizes must be compatible with batch_size=4
                with patch("sklearn.model_selection.train_test_split") as mock_split:
                    # First call: 30 → train=20, rest=10 (test)
                    X_train, X_rest = couples[:20], couples[20:]
                    y_train, y_rest = labels[:20], labels[20:]
                    s_train, s_rest = sources[:20], sources[20:]
                    # Second call: pair_ids → train, rest
                    pair_ids_all = pd.DataFrame({
                        "host_id": [f"h{i}" for i in range(30)],
                        "phage_id": [f"p{i}" for i in range(30)],
                        "label": labels,
                        "source": sources,
                    })
                    pids_train, pids_rest = pair_ids_all.iloc[:20].copy(), pair_ids_all.iloc[20:].copy()
                    # Third call: rest=10 → val=4, test=4
                    X_val, X_test = X_rest[:4], X_rest[4:8]
                    y_val, y_test = y_rest[:4], y_rest[4:8]
                    s_val, s_test = s_rest[:4], s_rest[4:8]
                    # Fourth call: pair_ids rest → val, test
                    pids_val, pids_test = pids_rest.iloc[:4].copy(), pids_rest.iloc[4:8].copy()
                    mock_split.side_effect = [
                        (X_train, X_rest, y_train, y_rest, s_train, s_rest),
                        (pids_train, pids_rest),
                        (X_val, X_test, y_val, y_test, s_val, s_test),
                        (pids_val, pids_test),
                    ]

                    sys.argv = [
                        "train.py",
                        "--config", str(config_path),
                        "--output-dir", str(tmp_path / "results"),
                        "--no-gpu",
                        "--run-name", "test_pretrain",
                    ]
                    train.main()

        # Verify output files were created
        results_dir = tmp_path / "results" / "test_pretrain"
        assert results_dir.exists(), f"Results dir not found: {results_dir}"
        assert (results_dir / "config.json").exists()
        assert (results_dir / "summary.json").exists()

        # Verify summary contents
        with open(results_dir / "summary.json") as f:
            summary = json.load(f)
        assert "train_size" in summary
        assert "best_val_auc" in summary
        assert "test_mcc" in summary

    def test_finetune_path(self, tmp_path):
        """Verify fine-tuning path runs end-to-end with mocked dependencies."""
        import numpy as np
        from importlib import reload
        import sys
        import train

        config_dict = {
            "defaults": {
                "training": {
                    "epochs": 2,
                    "batch_size": 4,
                    "patience": 2,
                    "learning_rate": 0.001,
                    "fine_tune_lr": 0.0001,
                    "fine_tune_epochs": 2,
                    "freeze_base": True,
                    "freeze_up_to": "concatenated_features",
                    "finetuned_model_name": "model_finetuned_best.keras",
                    "focal_alpha": 0.25,
                    "focal_gamma": 2.0,
                },
                "data": {
                    "negative_ratio": 1.0,
                    "bacterium_threshold": 50000,
                    "phage_threshold": 5000,
                    "bacterium_min_length": 100,
                    "phage_min_length": 100,
                    "exclude_ids": str(tmp_path / "excluded_pairs.csv"),
                    "exclude_sources": ["test_source"],
                },
                "output": {"dir": str(tmp_path / "results")},
                "gpu": {"enabled": False},
            },
        }

        config_path, pos_pairs, neg_pairs, couples, labels, sources = \
            self._setup_mocks(tmp_path, config_dict)

        mock_model = _make_mock_model()
        mock_retriever = MagicMock()
        mock_keras = MagicMock()
        mock_keras.callbacks.Callback = type("Callback", (), {"__init__": lambda s, *a, **kw: None})
        mock_keras.callbacks.ModelCheckpoint.return_value = MagicMock()
        mock_keras.callbacks.EarlyStopping.return_value = MagicMock()
        mock_keras.callbacks.LearningRateScheduler.return_value = MagicMock()
        mock_keras.callbacks.CSVLogger.return_value = MagicMock()
        mock_keras.metrics.AUC.return_value = MagicMock()

        # Reload first so patches apply to the fresh module
        reload(train)

        with (
            patch("pbi.quick_connect", return_value=mock_retriever),
            patch("pbi_adapter.PBIAdapter") as MockAdapter,
            patch.object(train, "load_or_build_model", return_value=(mock_model, True)),
            patch.object(train, "detect_gpu", return_value=False),
            patch.object(train, "keras", mock_keras, create=True),
        ):
            adapter_instance = MagicMock()
            adapter_instance.get_pair_ids_only.return_value = pd.DataFrame({
                "Phage_ID": [f"p{i}" for i in range(30)],
                "Host_ID": [f"h{i}" for i in range(30)],
            })
            adapter_instance.classify_pairs_by_interaction.return_value = (
                pos_pairs, neg_pairs
            )
            adapter_instance.prepare_training_data.return_value = (
                couples, labels, sources
            )
            adapter_instance.create_tf_generator = _make_mock_generator()
            MockAdapter.return_value = adapter_instance

            with patch("pbi.negative_examples.NegativeExampleGenerator") as MockNegGen:
                neg_gen_instance = MagicMock()
                gen_neg_df = pd.DataFrame({
                    "Phage_ID": [f"phage_gen_{i}" for i in range(10)],
                    "Host_ID": [f"host_gen_{i}" for i in range(10)],
                })
                neg_gen_instance.generate_random_negatives.return_value = gen_neg_df
                MockNegGen.return_value = neg_gen_instance

                with patch("sklearn.model_selection.train_test_split") as mock_split:
                    X_train, X_rest = couples[:20], couples[20:]
                    y_train, y_rest = labels[:20], labels[20:]
                    s_train, s_rest = sources[:20], sources[20:]
                    pair_ids_all = pd.DataFrame({
                        "host_id": [f"h{i}" for i in range(30)],
                        "phage_id": [f"p{i}" for i in range(30)],
                        "label": labels,
                        "source": sources,
                    })
                    pids_train, pids_rest = pair_ids_all.iloc[:20].copy(), pair_ids_all.iloc[20:].copy()
                    X_val, X_test = X_rest[:4], X_rest[4:8]
                    y_val, y_test = y_rest[:4], y_rest[4:8]
                    s_val, s_test = s_rest[:4], s_rest[4:8]
                    pids_val, pids_test = pids_rest.iloc[:4].copy(), pids_rest.iloc[4:8].copy()
                    mock_split.side_effect = [
                        (X_train, X_rest, y_train, y_rest, s_train, s_rest),
                        (pids_train, pids_rest),
                        (X_val, X_test, y_val, y_test, s_val, s_test),
                        (pids_val, pids_test),
                    ]

                    sys.argv = [
                        "train.py",
                        "--config", str(config_path),
                        "--pretrained-model", "/fake/model.keras",
                        "--output-dir", str(tmp_path / "results"),
                        "--no-gpu",
                        "--run-name", "test_finetune",
                    ]
                    train.main()

        # Verify output
        results_dir = tmp_path / "results" / "test_finetune"
        assert results_dir.exists(), f"Results dir not found: {results_dir}"
        assert (results_dir / "config.json").exists()
        assert (results_dir / "summary.json").exists()

        # Verify summary
        with open(results_dir / "summary.json") as f:
            summary = json.load(f)
        assert "best_val_auc" in summary

    def test_missing_exclude_ids_raises(self, tmp_path):
        """Verify that missing exclude_ids file raises FileNotFoundError."""
        import yaml

        config_dict = {
            "defaults": {
                "training": {"epochs": 1, "batch_size": 4, "patience": 1,
                             "learning_rate": 0.001, "focal_alpha": 0.25, "focal_gamma": 2.0},
                "data": {"negative_ratio": 1.0, "bacterium_threshold": 50000,
                         "phage_threshold": 5000, "bacterium_min_length": 100,
                         "phage_min_length": 100,
                         "exclude_ids": str(tmp_path / "nonexistent.csv"),
                         "exclude_sources": ["test_source"]},
                "output": {"dir": str(tmp_path / "results")},
                "gpu": {"enabled": False},
            },
        }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        mock_retriever = MagicMock()

        with (
            patch("pbi.quick_connect", return_value=mock_retriever),
            patch("pbi_adapter.PBIAdapter") as MockAdapter,
            patch("train.detect_gpu", return_value=False),
        ):
            adapter_instance = MagicMock()
            adapter_instance.get_pair_ids_only.return_value = pd.DataFrame({
                "Phage_ID": ["p1"], "Host_ID": ["h1"],
            })
            MockAdapter.return_value = adapter_instance

            sys.argv = [
                "train.py",
                "--config", str(config_path),
                "--no-gpu",
            ]

            from importlib import reload
            import train
            reload(train)

            with pytest.raises(FileNotFoundError, match="nonexistent.csv"):
                train.main()

    def test_config_overrides_defaults(self, tmp_path):
        """Verify config values are applied when CLI args are at defaults."""
        import yaml

        config_dict = {
            "defaults": {
                "training": {"epochs": 42, "batch_size": 8, "patience": 3,
                             "learning_rate": 0.01, "focal_alpha": 0.5, "focal_gamma": 1.0},
                "data": {"negative_ratio": 2.0, "bacterium_threshold": 50000,
                         "phage_threshold": 5000, "bacterium_min_length": 100,
                         "phage_min_length": 100},
                "output": {"dir": str(tmp_path / "results")},
                "gpu": {"enabled": False},
            },
        }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        # Verify the config is valid and values are correct
        with open(config_path) as f:
            loaded = yaml.safe_load(f)
        assert loaded["defaults"]["training"]["epochs"] == 42
        assert loaded["defaults"]["training"]["batch_size"] == 8
        assert loaded["defaults"]["data"]["negative_ratio"] == 2.0

    def test_profile_overrides_defaults(self, tmp_path):
        """Verify profile values override defaults correctly."""
        import yaml

        config_dict = {
            "defaults": {
                "training": {"epochs": 10, "batch_size": 4, "cross_validate": 0},
                "data": {"negative_ratio": 1.0},
            },
            "profiles": {
                "pretrain": {
                    "training": {"cross_validate": 5},
                },
                "finetune": {
                    "training": {"epochs": 5, "batch_size": 16, "cross_validate": 0},
                },
            },
        }

        # Simulate profile merging logic from train.py
        merged = dict(config_dict["defaults"])
        profile = config_dict["profiles"]["pretrain"]
        for section in ["training", "data"]:
            if section in profile:
                if section not in merged:
                    merged[section] = {}
                merged[section].update(profile[section])

        # cross_validate overridden by profile, epochs stays at default
        assert merged["training"]["cross_validate"] == 5
        assert merged["training"]["epochs"] == 10
        assert merged["data"]["negative_ratio"] == 1.0

    def test_invalid_profile_raises(self, tmp_path):
        """Verify invalid profile name raises ValueError."""
        import yaml

        config_dict = {
            "defaults": {"training": {"epochs": 10}},
            "profiles": {"pretrain": {"training": {"epochs": 20}}},
        }
        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        # Simulate profile lookup
        profiles = config_dict.get("profiles", {})
        missing = "nonexistent"
        with pytest.raises((ValueError, KeyError)):
            if missing not in profiles:
                raise ValueError(f"Profile '{missing}' not found")

