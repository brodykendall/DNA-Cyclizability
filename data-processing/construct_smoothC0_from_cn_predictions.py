import pandas as pd
import numpy as np
from numpy import array

# INPUTS:

# data_folder_path = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/"
data_folder_path = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/"

training_data_name = "tiling"
# training_data_name = "random"
# training_data_name = "chrv"
# training_data_name = "CN"

data_name_list = [
    # "yeast_chrI_1bpresolution_subsequence50",
    # "yeast_chrII_1bpresolution_subsequence50",
    # "yeast_chrIII_1bpresolution_subsequence50",
    # "yeast_chrIV_1bpresolution_subsequence50",
    # "yeast_chrV_1bpresolution_subsequence50",
    # "yeast_chrVI_1bpresolution_subsequence50",
    # "yeast_chrVII_1bpresolution_subsequence50",
    # "yeast_chrVIII_1bpresolution_subsequence50",
    # "yeast_chrIX_1bpresolution_subsequence50",
    # "yeast_chrX_1bpresolution_subsequence50",
    # "yeast_chrXI_1bpresolution_subsequence50",
    # "yeast_chrXII_1bpresolution_subsequence50",
    # "yeast_chrXIII_1bpresolution_subsequence50",
    # "yeast_chrXIV_1bpresolution_subsequence50",
    # "yeast_chrXV_1bpresolution_subsequence50",
    # "yeast_chrXVI_1bpresolution_subsequence50",
    # "tiling_full_reconstructed",
    # "expanded_random_library",
    "expanded_random_library2"
]

# data_file_type = ".txt"
data_file_type = ".csv"

# sequence_column_name = "Sequence"
sequence_column_name = "sequence"


# Original:
# for data_name in data_name_list:
#     df = pd.read_csv(f"{data_folder_path}{data_name}{data_file_type}")
    
#     cn_avg_pred_df = pd.read_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/ir_lstm_cn_{training_data_name}_{data_name}_predictions.csv")
#     pred_C26 = cn_avg_pred_df["n=26"]
#     pred_C29 = cn_avg_pred_df["n=29"]
#     pred_C31 = cn_avg_pred_df["n=31"]

#     pred_C26_rev_comp = cn_avg_pred_df["n=26_rev_comp"]
#     pred_C29_rev_comp = cn_avg_pred_df["n=29_rev_comp"]
#     pred_C31_rev_comp = cn_avg_pred_df["n=31_rev_comp"]

#     if data_name == "tiling_full_reconstructed":
#         # k9_exclusion_indices = df.loc[(df["region_pos"] == 0) | (df["region_pos"] == 1) | (df["region_pos"] == 2) | (df["region_pos"] == 3) | 
#         #                               (df["region_pos"] == 1001) | (df["region_pos"] == 1002) | (df["region_pos"] == 1003) | (df["region_pos"] == 1004)].index.values
#         k11_exclusion_indices = df.loc[(df["region_pos"] == 0) | (df["region_pos"] == 1) | (df["region_pos"] == 2) | (df["region_pos"] == 3) | (df["region_pos"] == 4) |
#                                     (df["region_pos"] == 1000) | (df["region_pos"] == 1001) | (df["region_pos"] == 1002) | (df["region_pos"] == 1003) | (df["region_pos"] == 1004)].index.values
#     else:
#         # k9_exclusion_indices = [0,1,2,3,df.shape[0]-4,df.shape[0]-3,df.shape[0]-2,df.shape[0]-1]
#         k11_exclusion_indices = [0,1,2,3,4,df.shape[0]-5,df.shape[0]-4,df.shape[0]-3,df.shape[0]-2,df.shape[0]-1]

#     k10_smooth_C26 = array([(pred_C26[i-5]+pred_C26[i+5])/20 + sum(pred_C26[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k10_smooth_C29 = array([(pred_C29[i-5]+pred_C29[i+5])/20 + sum(pred_C29[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k10_smooth_C31 = array([(pred_C31[i-5]+pred_C31[i+5])/20 + sum(pred_C31[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

#     k10_smooth_C26_rev_comp = array([(pred_C26_rev_comp[i-5]+pred_C26_rev_comp[i+5])/20 + sum(pred_C26_rev_comp[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k10_smooth_C29_rev_comp = array([(pred_C29_rev_comp[i-5]+pred_C29_rev_comp[i+5])/20 + sum(pred_C29_rev_comp[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k10_smooth_C31_rev_comp = array([(pred_C31_rev_comp[i-5]+pred_C31_rev_comp[i+5])/20 + sum(pred_C31_rev_comp[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

#     k11_smooth_C26 = array([np.mean(pred_C26[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k11_smooth_C29 = array([np.mean(pred_C29[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k11_smooth_C31 = array([np.mean(pred_C31[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

#     k11_smooth_C26_rev_comp = array([np.mean(pred_C26_rev_comp[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k11_smooth_C29_rev_comp = array([np.mean(pred_C29_rev_comp[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
#     k11_smooth_C31_rev_comp = array([np.mean(pred_C31_rev_comp[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

#     w_floor = 1/(10.4-10)
#     w_ceil = 1/(11-10.4)

#     smooth_C26 = array((k10_smooth_C26*w_floor + k11_smooth_C26*w_ceil)/(w_floor + w_ceil))
#     smooth_C29 = array((k10_smooth_C29*w_floor + k11_smooth_C29*w_ceil)/(w_floor + w_ceil))
#     smooth_C31 = array((k10_smooth_C31*w_floor + k11_smooth_C31*w_ceil)/(w_floor + w_ceil))

#     smooth_C26_rev_comp = array((k10_smooth_C26_rev_comp*w_floor + k11_smooth_C26_rev_comp*w_ceil)/(w_floor + w_ceil))
#     smooth_C29_rev_comp = array((k10_smooth_C29_rev_comp*w_floor + k11_smooth_C29_rev_comp*w_ceil)/(w_floor + w_ceil))
#     smooth_C31_rev_comp = array((k10_smooth_C31_rev_comp*w_floor + k11_smooth_C31_rev_comp*w_ceil)/(w_floor + w_ceil))

#     prelim_smoothC0 = (smooth_C26 + smooth_C31 + smooth_C26_rev_comp + smooth_C31_rev_comp + smooth_C29 + smooth_C29_rev_comp)/6

#     df["smooth_C26"] = smooth_C26
#     df["smooth_C29"] = smooth_C29
#     df["smooth_C31"] = smooth_C31
#     df["smooth_C26_rev_comp"] = smooth_C26_rev_comp
#     df["smooth_C29_rev_comp"] = smooth_C29_rev_comp
#     df["smooth_C31_rev_comp"] = smooth_C31_rev_comp

#     df["smoothC0"] = prelim_smoothC0

#     df.to_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/ir_lstm_cn_{training_data_name}_{data_name}_smoothC0_10_11.csv", index=False)

#     print(f"Finished /Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/ir_lstm_cn_{training_data_name}_{data_name}_smoothC0_10_11.csv")




# Updated (original scale, no rev comp):
# mean_C26 = -0.1702373294847514
# std_C26 = 0.5761561889426974
# mean_C29 = -0.1730189136679293
# std_C29 = 0.5889335849536195
# mean_C31 = -0.20523567238490675
# std_C31 = 0.6526103123263609
for data_name in data_name_list:
    df = pd.read_csv(f"{data_folder_path}{data_name}{data_file_type}")
    
    cn_avg_pred_df = pd.read_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/ir_lstm_cn_{training_data_name}_{data_name}_predictions.csv")
    pred_C26 = cn_avg_pred_df["n=26"]
    pred_C29 = cn_avg_pred_df["n=29"]
    pred_C31 = cn_avg_pred_df["n=31"]

    pred_C26_rev_comp = cn_avg_pred_df["n=26_rev_comp"]
    pred_C29_rev_comp = cn_avg_pred_df["n=29_rev_comp"]
    pred_C31_rev_comp = cn_avg_pred_df["n=31_rev_comp"]

    # pred_C26 = pred_C26*std_C26 + mean_C26
    # pred_C29 = pred_C29*std_C29 + mean_C29
    # pred_C31 = pred_C31*std_C31 + mean_C31

    if data_name == "tiling_full_reconstructed":
        k11_exclusion_indices = df.loc[(df["region_pos"] == 0) | (df["region_pos"] == 1) | (df["region_pos"] == 2) | (df["region_pos"] == 3) | (df["region_pos"] == 4) |
                                    (df["region_pos"] == 1000) | (df["region_pos"] == 1001) | (df["region_pos"] == 1002) | (df["region_pos"] == 1003) | (df["region_pos"] == 1004)].index.values
    else:
        k11_exclusion_indices = [0,1,2,3,4,df.shape[0]-5,df.shape[0]-4,df.shape[0]-3,df.shape[0]-2,df.shape[0]-1]

    k10_smooth_C26 = array([(pred_C26[i-5]+pred_C26[i+5])/20 + sum(pred_C26[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k10_smooth_C29 = array([(pred_C29[i-5]+pred_C29[i+5])/20 + sum(pred_C29[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k10_smooth_C31 = array([(pred_C31[i-5]+pred_C31[i+5])/20 + sum(pred_C31[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

    k10_smooth_C26_rev_comp = array([(pred_C26_rev_comp[i-5]+pred_C26_rev_comp[i+5])/20 + sum(pred_C26_rev_comp[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k10_smooth_C29_rev_comp = array([(pred_C29_rev_comp[i-5]+pred_C29_rev_comp[i+5])/20 + sum(pred_C29_rev_comp[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k10_smooth_C31_rev_comp = array([(pred_C31_rev_comp[i-5]+pred_C31_rev_comp[i+5])/20 + sum(pred_C31_rev_comp[(i-5+1):(i+5)])/10 if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

    k11_smooth_C26 = array([np.mean(pred_C26[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k11_smooth_C29 = array([np.mean(pred_C29[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k11_smooth_C31 = array([np.mean(pred_C31[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])

    k11_smooth_C26_rev_comp = array([np.mean(pred_C26_rev_comp[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k11_smooth_C29_rev_comp = array([np.mean(pred_C29_rev_comp[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])
    k11_smooth_C31_rev_comp = array([np.mean(pred_C31_rev_comp[(i-5):(i+5+1)]) if i not in k11_exclusion_indices else np.nan for i in range(df.shape[0])])


    w_floor = 1/(10.4-10)
    w_ceil = 1/(11-10.4)

    smooth_C26 = array((k10_smooth_C26*w_floor + k11_smooth_C26*w_ceil)/(w_floor + w_ceil))
    smooth_C29 = array((k10_smooth_C29*w_floor + k11_smooth_C29*w_ceil)/(w_floor + w_ceil))
    smooth_C31 = array((k10_smooth_C31*w_floor + k11_smooth_C31*w_ceil)/(w_floor + w_ceil))

    smooth_C26_rev_comp = array((k10_smooth_C26_rev_comp*w_floor + k11_smooth_C26_rev_comp*w_ceil)/(w_floor + w_ceil))
    smooth_C29_rev_comp = array((k10_smooth_C29_rev_comp*w_floor + k11_smooth_C29_rev_comp*w_ceil)/(w_floor + w_ceil))
    smooth_C31_rev_comp = array((k10_smooth_C31_rev_comp*w_floor + k11_smooth_C31_rev_comp*w_ceil)/(w_floor + w_ceil))

    prelim_smoothC0 = (smooth_C26 + smooth_C31 + smooth_C29 + smooth_C26_rev_comp + smooth_C29_rev_comp + smooth_C31_rev_comp)/6

    df["smooth_C26"] = smooth_C26
    df["smooth_C29"] = smooth_C29
    df["smooth_C31"] = smooth_C31
    df["smooth_C26_rev_comp"] = smooth_C26_rev_comp
    df["smooth_C29_rev_comp"] = smooth_C29_rev_comp
    df["smooth_C31_rev_comp"] = smooth_C31_rev_comp

    df["smoothC0"] = prelim_smoothC0

    df.to_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/ir_lstm_cn_{training_data_name}_{data_name}_smoothC0_10_11.csv", index=False)

    print(f"Finished /Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/ir_lstm_cn_{training_data_name}_{data_name}_smoothC0_10_11.csv")