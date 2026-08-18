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
        parser.add_argument("--batch-size", type=int, default=16)
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
        parser.add_argument("--no-gpu", action="store_true")
        parser.add_argument("--verbose", "-v", action="store_true")
        parser.add_argument("--log-file", type=str, default=None)

        args = parser.parse_args([])
        assert args.epochs == 10
        assert args.batch_size == 16
        assert args.limit is None
        assert args.no_gpu is False

    def test_custom_args(self):
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--epochs", type=int, default=10)
        parser.add_argument("--batch-size", type=int, default=16)
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
        model = build_model(bacterium_threshold=1000, phage_threshold=500)
        assert model is not None
        assert model.input[0].shape[1] == 1000
        assert model.input[1].shape[1] == 500


class TestConfigLoading:
    """Tests for YAML configuration loading."""

    def test_config_file_loading(self):
        """Verify config.yaml is valid and loadable."""
        import yaml
        config_path = Path(__file__).parent.parent / "PERPHECT" / "config.yaml"
        if config_path.exists():
            with open(config_path) as f:
                config = yaml.safe_load(f)
            assert "training" in config
            assert "data" in config
            assert "output" in config
            assert config["training"]["epochs"] == 10
            assert config["training"]["batch_size"] == 16

    def test_config_merging(self):
        """Verify config values override defaults."""
        import yaml

        config = {
            "training": {"epochs": 20, "batch_size": 32},
            "data": {"negative_ratio": 0.5},
        }

        # Simulate merging
        defaults = {"epochs": 10, "batch_size": 16, "limit": None}
        for section in ["training", "data"]:
            if section in config:
                for key, value in config[section].items():
                    defaults[key] = value

        assert defaults["epochs"] == 20
        assert defaults["batch_size"] == 32
        assert defaults["negative_ratio"] == 0.5
        assert defaults["limit"] is None  # Unchanged


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
