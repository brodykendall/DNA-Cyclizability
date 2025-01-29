import keras
import pandas as pd
import numpy as np
from numpy import array

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

# INPUTS:

# data_folder_path = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/"
data_folder_path = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/"

training_data_name = "tiling"
# training_data_name = "random"
# training_data_name = "chrv"
# training_data_name = "CN"

data_name_list = [
    # "ir_lstm_cn_tiling_yeast_chrI_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrII_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrIII_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrIV_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrVI_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrVII_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrVIII_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrIX_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrX_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrXI_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrXII_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrXIII_1bpresolution_subsequence50_smoothC0_10_11", 
    # "ir_lstm_cn_tiling_yeast_chrXIV_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrXV_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_yeast_chrXVI_1bpresolution_subsequence50_smoothC0_10_11",
    # "ir_lstm_cn_tiling_tiling_full_reconstructed_smoothC0_10_11",
    # "context_library",
    # "cycle1",
    # "cycle3",
    # "cycle5",
    # "cycle6",
    # "ir_lstm_cn_tiling_expanded_random_library_smoothC0_10_11", 
    "ir_lstm_cn_tiling_expanded_random_library2_smoothC0_10_11"
]



# data_name = "nuc_rotated" # csv
# data_name = "random_rotated" # csv
# data_name = "tiling_rotated" # csv
# data_name = "chrv_rotated" # csv

# data_file_type = ".txt"
data_file_type = ".csv"

# sequence_column_name = "Sequence"
sequence_column_name = "sequence"

training_column_name = "smoothC0"

model_identifier = "smoothC0_10_11_contracted"

# target_column_name_list = ["C0", "n=26", "n=29", "n=31"]
# target_column_name_list = ["C0", "C26", "C29", "C31", "smoothC0"]
target_column_name_list = ["smoothC0"]
# target_column_name_list = []

# Tiling contracted best fold - 9:
model = keras.models.load_model(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/benchmarks/ir_lstm_{model_identifier}_{training_data_name}_9.h5")

detrend_int = 0.001641373848542571
detrend_slope = 1.0158132314682007
# Mean and stdev of smoothC0 for contracted tiling library:
normal_mean = -0.011196041799376931
normal_std = 0.651684644408004

for data_name in data_name_list:
    df = pd.read_csv(f"{data_folder_path}{data_name}{data_file_type}")

    X = []
    for sequence_nt in df[sequence_column_name]:
        X.append(dnaOneHot(sequence_nt))
    X = array(X)
    X = X.reshape((X.shape[0],50,4,1))
    X_reverse = np.flip(X,[1,2])

    # Best fold:
    model_pred = model.predict(X)
    model_pred_reverse = model.predict(X_reverse)
    avg_predictions = (model_pred+model_pred_reverse)/2
    print(f"{training_column_name} prediction mean, std on {data_name} (pre-detrend): {np.mean(avg_predictions)}, {np.std(avg_predictions)}")
    avg_predictions = detrend_int + avg_predictions*detrend_slope
    avg_predictions = avg_predictions.flatten()
    print(f"{training_column_name} prediction mean, std on {data_name} (post-detrend): {np.mean(avg_predictions)}, {np.std(avg_predictions)}")

    for target_column_name in target_column_name_list:
        corr_df1 = pd.DataFrame({"target": df[target_column_name], "predicted": avg_predictions})
        # corr_df2 = pd.DataFrame({"target": df[target_column_name], "predicted": avg_predictions2})
        print(f"{training_column_name} prediction correlation on {data_name}, {target_column_name}: {corr_df1.corr().iloc[0,1]}")
        # print(f"{training_column_name} prediction correlation on {data_name}, {target_column_name}: {np.corrcoef(avg_predictions, df[target_column_name])[0,1]}"

    df[f"{training_column_name}_predictions"] = avg_predictions
    # df[f"{training_column_name}_predictions_unnorm"] = avg_predictions2

    df.to_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/predictions/ir_lstm_{model_identifier}_{training_data_name}_{data_name}_best_fold_predictions.csv", index=False)

    print(f"Finished /Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/predictions/ir_lstm_{model_identifier}_{training_data_name}_{data_name}_best_fold_predictions.csv\n\n")
