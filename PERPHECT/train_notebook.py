#!/usr/bin/env python3
"""
Auto-generated from PBI_Perphect_Training.ipynb

Run with: python {Path(py_path).name}
"""

# --- Cell 1 ---
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure PERPHECT directory is in path for local imports
sys.path.insert(0, os.path.dirname(os.path.abspath('__file__')))

from pbi import quick_connect
from pbi.negative_examples import NegativeExampleGenerator
from pbi_adapter import PBIAdapter
from transforms import translate_sequence_onehot

# --- Cell 2 ---
import tensorflow as tf

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    print(f'GPU available: {len(gpus)} device(s)')
    for gpu in gpus:
        print(f'  - {gpu.name}')
    print('Training will use GPU (fast).')
else:
    print('No GPU detected. Training will use CPU (slow).')
    print('See README.md for GPU setup instructions.')

# --- Cell 3 ---
retriever = quick_connect()
print(f"Connected to database.")
print(f"Has host data: {retriever._has_host_data}")

# --- Cell 4 ---
# Get summary statistics
phage_meta = retriever.get_phage_metadata()
print(f"Total phages: {len(phage_meta)}")
print(f"\nPhage length distribution:")
print(phage_meta['Length'].describe())

# --- Cell 5 ---
host_meta = retriever.get_host_metadata()
print(f"Total hosts: {len(host_meta)}")
print(f"\nHost genome length distribution:")
print(host_meta['Genome_Length'].describe())

# --- Cell 6 ---
# Load a small sample of positive pairs for initial testing
positive_pairs = retriever.get_phage_host_pairs(
    limit=1000,
    host_contig_mode='concat'
)
print(f"Loaded {len(positive_pairs)} positive pairs")
positive_pairs.head()

# --- Cell 7 ---
neg_gen = NegativeExampleGenerator(retriever)

# Generate equal number of negatives
negative_pairs = neg_gen.generate_random_negatives(
    positive_pairs,
    ratio=1.0
)
print(f"Generated {len(negative_pairs)} negative pairs")
negative_pairs.head()

# --- Cell 8 ---
adapter = PBIAdapter(
    retriever,
    bacterium_threshold=7_000_000,  # 7M bp (PERPHECT default)
    phage_threshold=200_000,        # 200K bp (PERPHECT default)
    bacterium_min_length=150_000,   # Filter smaller bacteria
    phage_min_length=1_500          # Filter smaller phages
)

couples_df, bacteria_df, phages_df = adapter.to_perphect_dataframes(
    positive_pairs,
    negative_pairs
)

print(f"Couples: {len(couples_df)}")
print(f"Bacteria: {len(bacteria_df)}")
print(f"Phages: {len(phages_df)}")
print(f"\nClass distribution:")
print(couples_df['interaction_type'].value_counts())

# --- Cell 9 ---
# Preview the data
print("Couples (first 5):")
display(couples_df.head())

print(f"\nBacteria sequences shape: {bacteria_df['bacterium_sequence'].iloc[0].shape}")
print(f"Phage sequences shape: {phages_df['phage_sequence'].iloc[0].shape}")

# --- Cell 10 ---
from sklearn.model_selection import train_test_split

# Prepare integer-indexed couples and labels
couples, labels = adapter.prepare_training_data(
    positive_pairs,
    negative_pairs
)

print(f"Total pairs: {len(couples)}")
print(f"Positive: {int(labels.sum())}, Negative: {int(len(labels) - labels.sum())}")

# Train/validation/test split (stratified)
X_train, X_test, y_train, y_test = train_test_split(
    couples, labels, stratify=labels, test_size=0.3, shuffle=True, random_state=42
)
X_valid, X_test, y_valid, y_test = train_test_split(
    X_test, y_test, stratify=y_test, test_size=0.5, shuffle=True, random_state=42
)

print(f"\nTrain: {len(X_train)}, Valid: {len(X_valid)}, Test: {len(X_test)}")

# --- Cell 11 ---
import math
import tensorflow as tf
import keras
from keras.models import Model
from keras.layers import Input, Conv1D, MaxPooling1D, Flatten, Dense, Dropout, Concatenate

BACTERIUM_THRESHOLD = 7_000_000
PHAGE_THRESHOLD = 200_000

# Bacteria branch
input1 = Input(shape=(BACTERIUM_THRESHOLD, 4), name="bacterial_input")
conv1_1 = Conv1D(64, 30, strides=10, activation='relu', name='bacterial_conv_1')(input1)
maxpool1_1 = MaxPooling1D(15, strides=5, name='bacterial_maxpool_1')(conv1_1)
conv1_2 = Conv1D(32, 25, strides=10, activation='relu', name='bacterial_conv_2')(maxpool1_1)
maxpool1_2 = MaxPooling1D(10, strides=5, name='bacterial_maxpool_2')(conv1_2)
conv1_3 = Conv1D(32, 10, strides=5, activation='relu', name='bacterial_conv_3')(maxpool1_2)
maxpool1_3 = MaxPooling1D(2, strides=2, name='bacterial_maxpool_3')(conv1_3)
flatten_bact = Flatten(name='bacteria_features')(maxpool1_3)

# Phage branch
input2 = Input(shape=(PHAGE_THRESHOLD, 4), name="phage_input")
conv2_1 = Conv1D(64, 30, strides=10, activation='relu', name='phage_conv_1')(input2)
maxpool2_1 = MaxPooling1D(15, strides=5, name='phage_maxpool_1')(conv2_1)
conv2_2 = Conv1D(32, 25, strides=10, activation='relu', name='phage_conv_2')(maxpool2_1)
maxpool2_2 = MaxPooling1D(2, strides=2, name='phage_maxpool_2')(conv2_2)
flatten_phage = Flatten(name='phage_features')(maxpool2_2)

# Classification head
concat_features = Concatenate(name='concatenated_features')([flatten_bact, flatten_phage])
dense1 = Dense(100, activation='relu')(concat_features)
dropout1 = Dropout(0.10)(dense1)
dense2 = Dense(1, activation='sigmoid')(dropout1)

model = Model(name='Perphect', inputs=[input1, input2], outputs=dense2)

# Learning rate schedule
def step_decay(epoch):
    initial_lrate = 0.0004
    drop = 0.5
    epochs_drop = 4.0
    return initial_lrate * math.pow(drop, math.floor((1 + epoch) / epochs_drop))

optimizer = keras.optimizers.Adam()
model.compile(optimizer, 'binary_crossentropy', metrics=['accuracy'])
model.summary()

# --- Cell 12 ---
from sklearn.utils import class_weight

BATCH_SIZE = 16
EPOCHS = 3
STEPS_PER_EPOCH = 400

# Create generators
train_gen = adapter.create_tf_generator(
    X_train, y_train, batch_size=BATCH_SIZE, shuffle=True
)
valid_gen = adapter.create_tf_generator(
    X_valid, y_valid, batch_size=BATCH_SIZE, shuffle=False
)
valid_steps = math.ceil(len(X_valid) / BATCH_SIZE)

# Callbacks
callbacks = [
    keras.callbacks.EarlyStopping(
        monitor='val_loss', patience=5, restore_best_weights=True
    ),
    keras.callbacks.LearningRateScheduler(step_decay),
]

# Train (class_weight not supported with generators in Keras 3.x;
# data is already balanced via ratio=1.0 in negative generation)
history = model.fit(
    train_gen,
    steps_per_epoch=STEPS_PER_EPOCH,
    epochs=EPOCHS,
    validation_data=valid_gen,
    validation_steps=valid_steps,
    callbacks=callbacks
)

# --- Cell 13 ---
from plotting_utils import historic_plots

acc_fig, loss_fig = historic_plots(history)
acc_fig.axes[0].set_title('Accuracy')
loss_fig.axes[0].set_title('Loss')

display(acc_fig)
display(loss_fig)

# Save figures
acc_fig.savefig('accuracy.png', dpi=150, bbox_inches='tight')
loss_fig.savefig('val_loss.png', dpi=150, bbox_inches='tight')
print('Saved accuracy.png and val_loss.png')

# --- Cell 14 ---
from sklearn.metrics import classification_report, confusion_matrix

# Test predictions
test_gen = adapter.create_tf_generator(
    X_test, y_test, batch_size=BATCH_SIZE, shuffle=False
)
test_steps = math.ceil(len(X_test) / BATCH_SIZE)
test_predictions = model.predict(test_gen, steps=test_steps)
test_pred_labels = (test_predictions.flatten() > 0.5).astype(int)

# Classification report
print('Classification Report:')
print(classification_report(y_test, test_pred_labels, target_names=['Negative', 'Positive']))

# Confusion matrix
cm = confusion_matrix(y_test, test_pred_labels)
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
ax.set(xticks=[0, 1], yticks=[0, 1],
       xticklabels=['Negative', 'Positive'],
       yticklabels=['Negative', 'Positive'],
       title='Confusion Matrix',
       ylabel='True label', xlabel='Predicted label')
plt.colorbar(im, ax=ax)
for i in range(2):
    for j in range(2):
        ax.text(j, i, format(cm[i, j], 'd'),
                ha='center', va='center',
                color='white' if cm[i, j] > cm.max() / 2 else 'black')
plt.tight_layout()
plt.show()

# --- Cell 15 ---
import os

output_dir = 'results'
os.makedirs(output_dir, exist_ok=True)

# Save model
model.save(os.path.join(output_dir, 'model.h5'))

# Save test results
results_df = pd.DataFrame({
    'bacterium_id': X_test[:, 0],
    'phage_id': X_test[:, 1],
    'observations': y_test,
    'predictions': test_predictions.flatten(),
    'prediction_labels': test_pred_labels
})
results_df.to_csv(os.path.join(output_dir, 'results_test_set.csv'), index=False)

print(f"Model saved to {output_dir}/model.h5")
print(f"Results saved to {output_dir}/results_test_set.csv")
