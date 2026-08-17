"""
Tests for PERPHECT plotting_utils module (historic_plots).
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).parent.parent / "PERPHECT"))


class TestHistoricPlots:
    """Tests for the historic_plots function."""

    def _make_history(self, metrics):
        """Create a mock Keras History object."""
        return SimpleNamespace(history=metrics)

    def test_accuracy_and_loss(self):
        from plotting_utils import historic_plots

        history = self._make_history({
            "accuracy": [0.6, 0.7, 0.8],
            "val_accuracy": [0.55, 0.65, 0.75],
            "loss": [0.5, 0.4, 0.3],
            "val_loss": [0.55, 0.45, 0.35],
        })

        acc_fig, loss_fig = historic_plots(history)

        assert acc_fig is not None
        assert loss_fig is not None
        assert len(acc_fig.axes) == 1
        assert len(loss_fig.axes) == 1

        # Check labels
        assert acc_fig.axes[0].get_xlabel() == "Epoch"
        assert acc_fig.axes[0].get_ylabel() == "Accuracy"
        assert loss_fig.axes[0].get_xlabel() == "Epoch"
        assert loss_fig.axes[0].get_ylabel() == "Loss"

        # Close figures to free memory
        import matplotlib.pyplot as plt
        plt.close(acc_fig)
        plt.close(loss_fig)

    def test_train_only(self):
        from plotting_utils import historic_plots

        history = self._make_history({
            "accuracy": [0.6, 0.7],
            "loss": [0.5, 0.4],
        })

        acc_fig, loss_fig = historic_plots(history)
        assert acc_fig is not None
        assert loss_fig is not None

        import matplotlib.pyplot as plt
        plt.close(acc_fig)
        plt.close(loss_fig)

    def test_empty_history(self):
        from plotting_utils import historic_plots

        history = self._make_history({})
        acc_fig, loss_fig = historic_plots(history)
        assert acc_fig is not None
        assert loss_fig is not None

        import matplotlib.pyplot as plt
        plt.close(acc_fig)
        plt.close(loss_fig)
