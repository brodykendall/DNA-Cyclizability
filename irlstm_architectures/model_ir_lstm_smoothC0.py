import time
import keras

from keras.models import Model
from keras.layers import Dropout, Flatten, BatchNormalization, TimeDistributed, Input, Add, Concatenate
from keras.layers import Dense, Conv2D, MaxPooling2D, LSTM, TimeDistributed, Reshape
import keras.backend as K
import keras.callbacks as callbacks

import pandas as pd
import numpy as np
from numpy import array
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

save_path = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/benchmarks/"
model_name = "ir_lstm_smoothC0_v1"
kf = KFold(n_splits = 10, shuffle =True)
num_epochs = 60

#### define functions ####

def model_cycle():
    inputs = Input(shape=(50, 4, 1))
        
    x = Conv2D(48, kernel_size=(3,4),
                   activation='relu',
                   padding='valid')(inputs)
    x = MaxPooling2D((2,1),padding='same')(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    # parallel line 1
    fx1 = Conv2D(48, kernel_size=(3,1),
                   activation='relu',
                   padding='same')(x)
    fx1 = BatchNormalization()(fx1)
    fx1 = Dropout(0.2)(fx1)
    fx1 = Conv2D(48, kernel_size=(3,1),
                   activation='relu',
                   padding='same')(fx1)
    fx1 = MaxPooling2D((2,1),padding='same')(fx1)
    fx1 = BatchNormalization()(fx1)
    fx1 = Dropout(0.2)(fx1)
    
    # parallel line 2
    fx2 = Conv2D(48, kernel_size=(11,1),
                   activation='relu',
                   padding='same')(x)
    fx2 = BatchNormalization()(fx2)
    fx2 = Dropout(0.2)(fx2)
    fx2 = Conv2D(48, kernel_size=(21,1),
                   activation='relu',
                   padding='same')(fx2)
    fx2 = MaxPooling2D((2,1),padding='same')(fx2)
    fx2 = BatchNormalization()(fx2)
    fx2 = Dropout(0.2)(fx2)
    
    # Add
    x1 = Concatenate(axis=-3)([fx1, fx2])
    x = Add()([x, x1])
    x = MaxPooling2D((2,1),padding='same')(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)
    
    x = Reshape((K.int_shape(x)[1], K.int_shape(x)[3]))(x)
    x = LSTM(20, return_sequences=False)(x)
    x = Dropout(0.2)(x)

    outputs = Dense(1, activation='linear')(x)
    network = Model(inputs, outputs)
    network.compile(optimizer='rmsprop',
                    loss='mean_squared_error')
    return network

def dnaOneHot(sequence):
    seq_array = array(list(sequence))
    code = {"A": [0], "C": [1], "G": [2], "T": [3], "N": [4],
            "a": [0], "c": [1], "g": [2], "t": [3], "n": [4]}
    onehot_encoded_seq = []
    for char in seq_array:
        onehot_encoded = np.zeros(5)
        onehot_encoded[code[char]] = 1
        onehot_encoded_seq.append(onehot_encoded[0:4])
    return onehot_encoded_seq

class TimeHistory(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.times = []

    def on_epoch_begin(self, batch, logs={}):
        self.epoch_time_start = time.process_time()

    def on_epoch_end(self, batch, logs={}):
        self.times.append(time.process_time() - self.epoch_time_start)

data_tiling = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/tiling_smoothC0.csv",delimiter = ",")
X5 = []
for sequence_nt in data_tiling["x50mer"]:
    X5.append(dnaOneHot(sequence_nt))
X5 = array(X5)
X5 = X5.reshape((X5.shape[0],50,4,1))
X5_reverse = np.flip(X5,[1,2])
Y5 = data_tiling["smooth_C0_v1"].values.astype(float)

data_chr5 = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/chrV_smoothC0.csv",delimiter = ",")
X6 = []
for sequence_nt in data_chr5["x50mer"]:
    X6.append(dnaOneHot(sequence_nt))
X6 = array(X6)
X6 = X6.reshape((X6.shape[0],50,4,1))
X6_reverse = np.flip(X6,[1,2])
Y6 = data_chr5["smooth_C0_v1"].values.astype(float)

m5 = np.mean(Y5)
std5 = np.std(Y5)
Z5 = (Y5-m5)/std5

m6 = np.mean(Y6)
std6 = np.std(Y6)
Z6 = (Y6-m6)/std6


#### tiling

VALIDATION_LOSS = []
fold_var = 1
n = Z5.shape[0]

fits = []
detrend = []
times = []
times2 = []

for train_index, val_index in kf.split(Z5):
    training_X = X5[train_index]
    training_X_reverse = X5_reverse[train_index]
    validation_X = X5[val_index]
    validation_X_reverse = X5_reverse[val_index]
    training_Y = Z5[train_index]
    validation_Y = Z5[val_index]
    # CREATE NEW MODEL
    model = model_cycle()
    # CREATE CALLBACKS
    checkpoint = callbacks.ModelCheckpoint(save_path + model_name+"_tiling_"+str(fold_var)+".h5",
                                                    monitor='val_loss', verbose=1,
                                                    save_best_only=True, mode='min')
    time_callback = TimeHistory()

    history = model.fit(training_X, training_Y,
                        epochs=num_epochs,
                        callbacks= [checkpoint, time_callback],
                        validation_data=(validation_X, validation_Y))
    model.load_weights(save_path + model_name+"_tiling_"+str(fold_var)+".h5")
    model.save(save_path+model_name+"_tiling_"+str(fold_var),save_traces=False)
    times.append(time_callback.times)

    pred_Y = model.predict(training_X)
    pred_Y = pred_Y.reshape(pred_Y.shape[0])
    pred_Y_reverse = model.predict(training_X_reverse)
    pred_Y_reverse = pred_Y_reverse.reshape(pred_Y_reverse.shape[0])
    pred_Y = (pred_Y+pred_Y_reverse)/2
    reg =  LinearRegression().fit(array(pred_Y).reshape(-1, 1), array(training_Y).reshape(-1, 1))
    
    detrend_int = reg.intercept_
    detrend_slope = reg.coef_
    detrend.append([float(detrend_int), float(detrend_slope)])
    
    start_time = time.process_time()
    fit = model.predict(validation_X)
    fit = fit.reshape(fit.shape[0])
    fit_reverse = model.predict(validation_X_reverse)
    fit_reverse = fit_reverse.reshape(fit_reverse.shape[0])
    reverse_corr = np.corrcoef(fit, fit_reverse)[0,1]
    fit = (fit + fit_reverse)/2
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, validation_Y)[0,1],np.mean(np.square(fit-validation_Y)),np.mean(fit),np.std(fit),reverse_corr]
    fits.append(fit_tmp)
    fit = detrend_int + fit * detrend_slope
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, validation_Y)[0,1],np.mean(np.square(fit-validation_Y)),np.mean(fit),np.std(fit),reverse_corr]
    time0 = time.process_time() - start_time
    times2.append([time0])
    fits.append(fit_tmp)
    
    start_time = time.process_time()
    fit = model.predict(X6)
    fit = fit.reshape(fit.shape[0])
    fit_reverse = model.predict(X6_reverse)
    fit_reverse = fit_reverse.reshape(fit_reverse.shape[0])
    reverse_corr = np.corrcoef(fit, fit_reverse)[0,1]
    fit = (fit + fit_reverse)/2
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, Z6)[0,1],np.mean(np.square(fit-Z6)),np.mean(fit),np.std(fit),reverse_corr]
    fits.append(fit_tmp)
    fit = detrend_int + fit * detrend_slope
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, Z6)[0,1],np.mean(np.square(fit-Z6)),np.mean(fit),np.std(fit),reverse_corr]
    time0 = time.process_time() - start_time
    times2.append([time0])
    fits.append(fit_tmp)
    
    K.clear_session()
    fold_var += 1
    
detrend = array(detrend)
detrend = pd.DataFrame(detrend)
detrend.to_csv(save_path +model_name+"_detrend_tiling.txt", index = False)

fits = array(fits)
fits = pd.DataFrame((fits))
fits.to_csv(save_path +model_name+ "_fits_tiling.txt", index = False)

with open(save_path +model_name+"_time_tiling.txt", "w") as file:
    for row in times:
        s = " ".join(map(str, row))
        file.write(s+'\n')

with open(save_path +model_name+"_pred_time_tiling.txt", "w") as file:
    for row in times2:
        s = " ".join(map(str, row))
        file.write(s+'\n')





#### chrv

VALIDATION_LOSS = []
fold_var = 1
n = Z6.shape[0]

fits = []
detrend = []
times = []
times2 = []

for train_index, val_index in kf.split(Z6):
    training_X = X6[train_index]
    training_X_reverse = X6_reverse[train_index]
    validation_X = X6[val_index]
    validation_X_reverse = X6_reverse[val_index]
    training_Y = Z6[train_index]
    validation_Y = Z6[val_index]
    # CREATE NEW MODEL
    model = model_cycle()
    # CREATE CALLBACKS
    checkpoint = callbacks.ModelCheckpoint(save_path + model_name+"_chrv_"+str(fold_var)+".h5",
                                                    monitor='val_loss', verbose=1,
                                                    save_best_only=True, mode='min')
    time_callback = TimeHistory()

    history = model.fit(training_X, training_Y,
                        epochs=num_epochs,
                        callbacks= [checkpoint, time_callback],
                        validation_data=(validation_X, validation_Y))
    model.load_weights(save_path + model_name+"_chrv_"+str(fold_var)+".h5")
    model.save(save_path+ model_name+"_chrv_"+str(fold_var),save_traces=False)
    times.append(time_callback.times)

    pred_Y = model.predict(training_X)
    pred_Y = pred_Y.reshape(pred_Y.shape[0])
    pred_Y_reverse = model.predict(training_X_reverse)
    pred_Y_reverse = pred_Y_reverse.reshape(pred_Y_reverse.shape[0])
    pred_Y = (pred_Y+pred_Y_reverse)/2
    reg =  LinearRegression().fit(array(pred_Y).reshape(-1, 1), array(training_Y).reshape(-1, 1))
    
    detrend_int = reg.intercept_
    detrend_slope = reg.coef_
    detrend.append([float(detrend_int), float(detrend_slope)])
    
    start_time = time.process_time()
    fit = model.predict(X5)
    fit = fit.reshape(fit.shape[0])
    fit_reverse = model.predict(X5_reverse)
    fit_reverse = fit_reverse.reshape(fit_reverse.shape[0])
    reverse_corr = np.corrcoef(fit, fit_reverse)[0,1]
    fit = (fit + fit_reverse)/2
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, Z5)[0,1],np.mean(np.square(fit-Z5)),np.mean(fit),np.std(fit),reverse_corr]
    fits.append(fit_tmp)
    fit = detrend_int + fit * detrend_slope
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, Z5)[0,1],np.mean(np.square(fit-Z5)),np.mean(fit),np.std(fit),reverse_corr]
    time0 = time.process_time() - start_time
    times2.append([time0])
    fits.append(fit_tmp)
    
    start_time = time.process_time()
    fit = model.predict(validation_X)
    fit = fit.reshape(fit.shape[0])
    fit_reverse = model.predict(validation_X_reverse)
    fit_reverse = fit_reverse.reshape(fit_reverse.shape[0])
    reverse_corr = np.corrcoef(fit, fit_reverse)[0,1]
    fit = (fit + fit_reverse)/2
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, validation_Y)[0,1],np.mean(np.square(fit-validation_Y)),np.mean(fit),np.std(fit),reverse_corr]
    fits.append(fit_tmp)
    fit = detrend_int + fit * detrend_slope
    fit = fit.flatten()
    fit_tmp =[np.corrcoef(fit, validation_Y)[0,1],np.mean(np.square(fit-validation_Y)),np.mean(fit),np.std(fit),reverse_corr]
    time0 = time.process_time() - start_time
    times2.append([time0])
    fits.append(fit_tmp)
    
    K.clear_session()
    fold_var += 1
    
detrend = array(detrend)
detrend = pd.DataFrame(detrend)
detrend.to_csv(save_path +model_name+"_detrend_chrv.txt", index = False)

fits = array(fits)
fits = pd.DataFrame((fits))
fits.to_csv(save_path +model_name+"_fits_chrv.txt", index = False)

with open(save_path +model_name+"_time_chrv.txt", "w") as file:
    for row in times:
        s = " ".join(map(str, row))
        file.write(s+'\n')

with open(save_path +model_name+"_pred_time_chrv.txt", "w") as file:
    for row in times2:
        s = " ".join(map(str, row))
        file.write(s+'\n')