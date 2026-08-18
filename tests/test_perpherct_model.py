"""
Smoke test for PERPHECT model build (Keras 3.x compatibility).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "PERPHECT"))


@pytest.fixture(scope="module")
def perphect_model():
    """Build the PERPHECT model and return it."""
    keras = pytest.importorskip("keras")
    from keras.models import Model
    from keras.layers import Input, Conv1D, MaxPooling1D, Flatten, Dense, Dropout, Concatenate
    import math

    BACTERIUM_THRESHOLD = 1000  # Minimum for simplified 2-layer architecture
    PHAGE_THRESHOLD = 500

    # Bacteria branch
    input1 = Input(shape=(BACTERIUM_THRESHOLD, 4), name="bacterial_input")
    conv1_1 = Conv1D(64, 30, strides=10, activation='relu', name='bacterial_conv_1')(input1)
    maxpool1_1 = MaxPooling1D(15, strides=5, name='bacterial_maxpool_1')(conv1_1)
    conv1_2 = Conv1D(32, 25, strides=10, activation='relu', name='bacterial_conv_2')(maxpool1_1)
    flatten_bact = Flatten(name='bacteria_features')(conv1_2)

    # Phage branch
    input2 = Input(shape=(PHAGE_THRESHOLD, 4), name="phage_input")
    conv2_1 = Conv1D(64, 30, strides=10, activation='relu', name='phage_conv_1')(input2)
    maxpool2_1 = MaxPooling1D(15, strides=5, name='phage_maxpool_1')(conv2_1)
    flatten_phage = Flatten(name='phage_features')(maxpool2_1)

    # Classification head
    concat_features = Concatenate(name='concatenated_features')([flatten_bact, flatten_phage])
    dense1 = Dense(100, activation='relu')(concat_features)
    dropout1 = Dropout(0.10)(dense1)
    dense2 = Dense(1, activation='sigmoid')(dropout1)

    model = Model(name='Perphect', inputs=[input1, input2], outputs=dense2)
    return model


class TestPERPHECTModel:
    """Tests for the PERPHECT Keras model."""

    def test_model_builds(self, perphect_model):
        """Model should build without errors."""
        assert perphect_model is not None

    def test_model_summary(self, perphect_model, capsys):
        """Model summary should print without errors."""
        perphect_model.summary()
        captured = capsys.readouterr()
        assert "Perphect" in captured.out

    def test_model_compile(self, perphect_model):
        """Model should compile without errors."""
        import keras
        perphect_model.compile(
            optimizer=keras.optimizers.Adam(),
            loss='binary_crossentropy',
            metrics=['accuracy']
        )

    def test_model_predict(self, perphect_model):
        """Model should predict without errors."""
        import numpy as np
        import keras

        perphect_model.compile(
            optimizer=keras.optimizers.Adam(),
            loss='binary_crossentropy',
            metrics=['accuracy']
        )

        # Small batch
        bacterium_input = np.random.rand(2, 1000, 4).astype(np.float32)
        phage_input = np.random.rand(2, 500, 4).astype(np.float32)

        predictions = perphect_model.predict(
            [bacterium_input, phage_input],
            verbose=0
        )

        assert predictions.shape == (2, 1)
        assert all(0 <= p <= 1 for p in predictions.flatten())
