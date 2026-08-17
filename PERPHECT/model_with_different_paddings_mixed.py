from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import os
import math

import tensorflow as tf
import keras
from keras.models import Model 
from keras import layers
from keras import Input
from plotting_utils import historic_plots
from sklearn.utils import class_weight
from transforms import translate_sequence_onehot


import argparse

BACTERIUM_THRESHOLD = 7000000
PHAGE_THRESHOLD = 200000

BACTERIUM_MIN_LENGTH = 150000
PHAGE_MIN_LENGTH = 1500

def main(batch_size, epochs, steps_per_epoch, train_filepath, output_filepath, zero_padded):

    if zero_padded:
        bacteria_filename = 'zero_padded_bacteria_df.csv'
        phages_filename = 'zero_padded_phages_df.csv'
    else:
        bacteria_filename = 'replicate_padded_bacteria_df.csv'
        phages_filename = 'replicate_padded_phages_df.csv'

    print('bacteria_filename = ', bacteria_filename)
    print('phages_filename = ', phages_filename)

    ####################################################################################################################

    print("Keras version : {}".format(keras.__version__))
    print("TensorFlow version: {}".format(tf.version.VERSION))

    # Keras callbacks
    callbacks_list = [

        # Save best model
        keras.callbacks.ModelCheckpoint(filepath=os.path.join(output_filepath, 'model.h5'),
                        monitor='val_loss', save_best_only=True),

        # Stop if validation loss does not decrease in several epochs (restore best though)
        keras.callbacks.EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True),

    ]
    
    ####################################################################################################################

    print("Read data...")

    couples_df = pd.read_csv(os.path.join(train_filepath, 'cleaned_couples_df.csv'))
    bacteria_df = pd.read_csv(os.path.join(train_filepath, bacteria_filename), index_col='bacterium_id')
    phages_df = pd.read_csv(os.path.join(train_filepath, phages_filename), index_col='phage_id')

    bacteria_df['bacterium_sequence_onehot'] = bacteria_df['bacterium_sequence'].astype(str).apply(translate_sequence_onehot)
    phages_df['phage_sequence_onehot'] = phages_df['phage_sequence'].astype(str).apply(translate_sequence_onehot)

    X = couples_df.loc[:,'bacterium_id':'phage_id'].values
    y = couples_df['interaction_type'].values

    print(couples_df.head())
    print('couples_df[interaction_type].unique() = ', couples_df['interaction_type'].unique())
    print(couples_df.groupby('interaction_type').size())

    # Split data into train and test sets
    # X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.5, shuffle=True, random_state=0)
    # X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, stratify=y_train, test_size=0.2, shuffle=True, random_state=0)


    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.7, shuffle=True, random_state=0)
    X_valid, X_test, y_valid, y_test = train_test_split(X_test, y_test, stratify=y_test, test_size=0.5, shuffle=True, random_state=0)

    print('x_train.shape = ', X_train.shape)
    print('X_valid.shape = ', X_valid.shape)
    print('X_test.shape = ', X_test.shape)
    print('y_train.shape = ', y_train.shape)
    print('y_valid.shape = ', y_valid.shape)
    print('y_test.shape = ', y_test.shape)

    ####################################################################################################################

    # Transform data
    print("Transform data...")

    print("\tTranslate bacteria...")
    bacterium_sequences = bacteria_df.to_dict().get('bacterium_sequence_onehot')

    print("\tTranslate phages...")
    phage_sequences = phages_df.to_dict().get('phage_sequence_onehot')

    print('len(bacterium_sequences) = ', len(bacterium_sequences))
    print('len(phage_sequences) = ', len(phage_sequences))

    ####################################################################################################################

    print("Building model...")
    class_weights = class_weight.compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)


    with keras.backend.name_scope('bacterial_convolution'):

        input1 = Input(shape=(BACTERIUM_THRESHOLD, 4), name="bacterial_input")

        # First convolution layer for input 1
        conv1_1 = layers.Conv1D(name="bacterial_conv_1",
                         # input_shape=(BACTERIA_THRESHOLD, 4),
                         activation="relu",
                         filters=64,
                         kernel_size=30,
                         strides=10,
                         # padding="causal"
                         # kernel_initializer=keras.initializers.Constant(value=0.25)
                         )(input1)

        maxpooling1_1 = layers.MaxPooling1D(name="bacterial_maxpool_1",
                                     pool_size=15,
                                     strides=5
                                     )(conv1_1)

        # Second convolution layer for input 1
        conv1_2 = layers.Conv1D(name="bacterial_conv_2",
                         activation="relu",
                         filters=32,
                         kernel_size=25,
                         strides=10,
                         # padding="causal"
                         )(maxpooling1_1)

        maxpooling1_2 = layers.MaxPooling1D(name="bacterial_maxpool_2",
                                     pool_size=10,
                                     strides=5)(conv1_2)

        # Third convolution layer for input 1
        conv1_3 = layers.Conv1D(name="bacterial_conv_3",
                         activation="relu",
                         filters=32,
                         kernel_size=10,
                         strides=5,
                         # padding='causal'
                         )(maxpooling1_2)

        maxpooling1_3 = layers.MaxPooling1D(name="bacterial_maxpool_3",
                                     pool_size=2,
                                     strides=2
                                     )(conv1_3)

        flatten_bact = layers.Flatten(name="bacteria_features")(maxpooling1_3)

    with keras.backend.name_scope('phage_convolution'):

        input2 = Input(shape=(PHAGE_THRESHOLD, 4), name="phage_input")

        # First convolution layer for input 1
        conv2_1 = layers.Conv1D(name="phage_conv_1",
                         # input_shape=(PHAGE_THRESHOLD, 4),
                         activation="relu",
                         filters=64,
                         kernel_size=30,
                         strides=10,
                         # padding="causal"
                         # kernel_initializer=keras.initializers.Constant(value=0.25)
                         )(input2)

        maxpooling2_1 = layers.MaxPooling1D(name="phage_maxpool_1",
                                     pool_size=15,
                                     strides=5
                                     )(conv2_1)

        # Second convolution layer for input 1
        conv2_2 = layers.Conv1D(name="phage_conv_2",
                         activation="relu",
                         filters=32,
                         kernel_size=25,
                         strides=10
                         )(maxpooling2_1)

        maxpooling2_2 = layers.MaxPooling1D(name="phage_maxpool_2",
                                     pool_size=2,
                                     strides=2
                                     )(conv2_2)

        flatten_phage = layers.Flatten(name="phage_features")(maxpooling2_2)

    concat_features = layers.Concatenate(name="concatenated_features")([flatten_bact, flatten_phage])

    # dropout1 = Dropout(0.10)(concat_features)
    dense1 = layers.Dense(100, activation="relu")(concat_features)
    dropout1 = layers.Dropout(0.10)(dense1)
    dense2 = layers.Dense(1, activation="sigmoid")(dropout1)

    model = Model(name="Perphect",
                  inputs=[input1, input2],
                  outputs=dense2)

    # Configure optimizer
    # learning_rate = 0.0001

    optimizer = keras.optimizers.Adam()

    # Learning rate step decay
    def step_decay(epoch):
        initial_lrate = 0.0004
        drop = 0.5
        epochs_drop = 4.0
        lrate = initial_lrate * math.pow(drop, math.floor((1 + epoch) / epochs_drop))
        print("lr", lrate)
        return lrate

    callbacks_list.append(keras.callbacks.LearningRateScheduler(step_decay))

    # Configuring and compiling model
    model.compile(optimizer, 'binary_crossentropy', metrics=['accuracy'])

    model.summary()

    ####################################################################################################################

    print("Training model...")

    # Model generators
    train_generator = generator(X_train, y_train, 
                                bacterium_sequences, phage_sequences, 
                                min_index=0, 
                                max_index=None, 
                                shuffle=True, 
                                batch_size=batch_size)

    valid_generator = generator(X_valid, y_valid, 
                                bacterium_sequences, phage_sequences, 
                                min_index=0, 
                                max_index=None, 
                                shuffle=False, 
                                batch_size=batch_size)

    # Steps per generator
    valid_steps = math.ceil(len(X_valid) / batch_size)

    hist = model.fit(train_generator,
                     steps_per_epoch=steps_per_epoch,
                     epochs=epochs,
                     class_weight=dict(zip(np.unique(y_train), class_weights)),
                     validation_data=valid_generator,
                     validation_steps=valid_steps,
                     callbacks=callbacks_list)

    # Plotting accuracy and loss
    acc_plt, loss_plt = historic_plots(hist)

    acc_plt.axes[0].set_title("Accuracy")
    acc_plt.savefig(os.path.join(output_filepath, 'accuracy.png'))

    loss_plt.axes[0].set_title("Loss")
    loss_plt.savefig(os.path.join(output_filepath, 'val_loss.png'))

    # Validation generator
    valid_generator = generator(X_valid, y_valid, 
                                bacterium_sequences, phage_sequences, 
                                min_index=0, 
                                max_index=None, 
                                shuffle=False, 
                                batch_size=batch_size)

    # Validation predict
    val_predictions = model.predict(valid_generator, steps=valid_steps)

    # Save results to a CSV file
    results_df = pd.DataFrame()
    results_df['bacterium_id'] = X_valid[:,0]
    results_df['phage_id'] = X_valid[:,1]
    results_df['observations'] = y_valid
    results_df['predictions'] = val_predictions
    results_df['prediction_labels'] = results_df['predictions'].apply(round)

    results_df.to_csv(os.path.join(output_filepath, 'results_validation_set.csv'), index=False)

    ####################################################################################################################

    print("Testing in vivo data...")

    X = X_test
    y = y_test

    # Generator of one_hot encoded padded sequences
    test_generator = generator(X, y, 
                                bacterium_sequences, phage_sequences, 
                                min_index=0, 
                                max_index=None, 
                                shuffle=False, 
                                batch_size=batch_size)

    # Steps per generator
    test_steps = math.ceil(len(X) / batch_size)

    # Test predict
    test_predictions = model.predict(test_generator, steps=test_steps)

    # Save results to a CSV file
    results_df = pd.DataFrame()
    results_df['bacterium_id'] = X[:,0]
    results_df['phage_id'] = X[:,1]
    results_df['observations'] = y
    results_df['predictions'] = test_predictions
    results_df['prediction_labels'] = results_df['predictions'].apply(round)

    results_df.to_csv(os.path.join(output_filepath, 'results_test_set.csv'), index=False)

def generator(couples, labels, bacterium_sequences, phage_sequences, min_index, max_index, shuffle=False, batch_size=16):
    if max_index is None:
        max_index = len(couples)
    i = min_index
    
    while True:
        if shuffle:
            rows = np.random.randint(min_index , max_index, size=batch_size)
        else:
            rows = np.arange(i, min(i + batch_size, max_index))
            i += len(rows)
            if(i >= max_index):
                i = min_index
    
        bacterium_samples = np.zeros((len(rows), BACTERIUM_THRESHOLD, 4), dtype=np.uint8)
        phage_samples = np.zeros((len(rows), PHAGE_THRESHOLD, 4), dtype=np.uint8)
        targets = np.zeros((len(rows),))
        
        for j, row in enumerate(rows):
            bacterium_samples[j] = bacterium_sequences[couples[row, 0]]
            phage_samples[j] = phage_sequences[couples[row, 1]]
            targets[j] = labels[row]
            
        yield [bacterium_samples, phage_samples], targets

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", help="Batch size", type=int, default=16)
    parser.add_argument("--epoch", help="Num epochs", type=int, default=3)

    steps_per_epoch = 400
    train_filepath = 'data/mixed_data_set_'
    output_filepath = 'results/mixed_train_set_zero_balanced_3'
    zero_padded = True

    args = parser.parse_args()

    main(args.batch, args.epoch, steps_per_epoch, train_filepath, output_filepath, zero_padded)

