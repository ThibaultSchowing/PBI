"""
Plotting utilities for PERPHECT model training visualization.

Provides functions to generate training history plots (accuracy, loss).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def historic_plots(history):
    """
    Generate accuracy and loss plots from Keras training history.

    Args:
        history: Keras History object (result of model.fit())

    Returns:
        Tuple of (accuracy_figure, loss_figure).
    """
    acc_fig, acc_ax = plt.subplots(1, 1, figsize=(8, 5))
    loss_fig, loss_ax = plt.subplots(1, 1, figsize=(8, 5))

    # Accuracy
    if 'accuracy' in history.history:
        acc_ax.plot(history.history['accuracy'], label='Train', linewidth=2)
    if 'val_accuracy' in history.history:
        acc_ax.plot(history.history['val_accuracy'], label='Validation', linewidth=2)
    acc_ax.set_xlabel('Epoch')
    acc_ax.set_ylabel('Accuracy')
    acc_ax.legend()
    acc_ax.grid(True, alpha=0.3)

    # Loss
    if 'loss' in history.history:
        loss_ax.plot(history.history['loss'], label='Train', linewidth=2)
    if 'val_loss' in history.history:
        loss_ax.plot(history.history['val_loss'], label='Validation', linewidth=2)
    loss_ax.set_xlabel('Epoch')
    loss_ax.set_ylabel('Loss')
    loss_ax.legend()
    loss_ax.grid(True, alpha=0.3)

    return acc_fig, loss_fig
