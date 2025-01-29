from train import *
from datasets import *
import sys
from sklearn.model_selection import KFold


hyperparameters = {"nodes_per_layer":16,"num_epochs":200, "learning_rate":0.00001, "batch_size":32,
                   "num_layers":5, "max_time":8*60*60, "dropout": 0.0, "patience": 30}

target_variable = "smooth_10.4bp_cn_mean2"

# yeast_chrV_matched = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/yeast_chrV_ir_lstm_tiling_post_smoothed_matched.csv")

data_tiling = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/tiling_ir_lstm_cn_tiling_post_smoothed.csv")
data_tiling = data_tiling.loc[~data_tiling["C26"].isna()]

print(data_tiling.shape[0])

X = np.array(data_tiling[["C26", "C29", "C31"]])
# X = np.array(data_tiling[["n=26", "n=29", "n=31"]])

# X = np.array(yeast_chrV_matched[["n=26", "n=29", "n=31"]])
# X = np.array(yeast_chrV_matched[["n=26", "n=29", "n=31", "ps1_fourier_amp", "ps1_fourier_phase"]])
# X = np.array(yeast_chrV_matched[["n=26", "n=29", "n=31", 
#                                  "A_fourier_abs_5", "A_fourier_abs_10", "A_fourier_phase_5", "A_fourier_phase_10",
#                                  "C_fourier_abs_5", "C_fourier_abs_10", "C_fourier_phase_5", "C_fourier_phase_10", 
#                                  "G_fourier_abs_5", "G_fourier_abs_10", "G_fourier_phase_5", "G_fourier_phase_10", 
#                                  "T_fourier_abs_5", "T_fourier_abs_10", "T_fourier_phase_5", "T_fourier_phase_10", 
#                                  "AT_fourier_abs_5", "AT_fourier_abs_10", "AT_fourier_phase_5", "AT_fourier_phase_10"]])
# X = np.array(yeast_chrV_matched[["n=26", "n=29", "n=31", 
#                                  "A_fourier_abs_2", "A_fourier_abs_5", "A_fourier_abs_10", 
#                                  "A_fourier_phase_2", "A_fourier_phase_5", "A_fourier_phase_10",
#                                  "C_fourier_abs_2", "C_fourier_abs_5", "C_fourier_abs_10", 
#                                  "C_fourier_phase_2", "C_fourier_phase_5", "C_fourier_phase_10", 
#                                  "G_fourier_abs_2", "G_fourier_abs_5", "G_fourier_abs_10", 
#                                  "G_fourier_phase_2", "G_fourier_phase_5", "G_fourier_phase_10", 
#                                  "T_fourier_abs_2", "T_fourier_abs_5", "T_fourier_abs_10", 
#                                  "T_fourier_phase_2", "T_fourier_phase_5", "T_fourier_phase_10", 
#                                  "AT_fourier_abs_2", "AT_fourier_abs_5", "AT_fourier_abs_10", 
#                                  "AT_fourier_phase_2", "AT_fourier_phase_5", "AT_fourier_phase_10"]])
# y = np.array(yeast_chrV_matched[["DNAcycP_pred_chrV_post_smooth_5bp"]])

y = np.array(data_tiling[[target_variable]])


X = X[~np.isnan(y).flatten()]
y = y[~np.isnan(y).flatten()]

logging_info = {"print_every":1, "folder_path":"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/predict_post_smoothed/",
                "file_name": f"5bp_DNAcycP_{hyperparameters['nodes_per_layer']}nodes_per_layer_{hyperparameters['num_layers']}layers_{int(hyperparameters['dropout']*10)}dropout_{X.shape[1]}input_features"}
log_file_path = logging_info["folder_path"] + "log_files/" + logging_info["file_name"] + ".txt"

with open(log_file_path, "w") as f:
    sys.stdout = f

    print(f"Hyperparameters: \n \t {hyperparameters}")
    
    kf = KFold(10, shuffle=True)

    fold=0
    for train_index, val_index in kf.split(X):

        print(f"\n \n \n Begin Training Fold {fold}: \n \n \n")

        training_data = PostSmoothedDataset(X[train_index], y[train_index])
        validation_data = PostSmoothedDataset(X[val_index], y[val_index])

        model, train_loss, val_loss, val_corr = train(training_data, hyperparameters=hyperparameters, 
                                                      logging_info=logging_info, model_name="FeedForward", val_data=validation_data, 
                                                      plot_0_val=False, fold=fold, in_features=X.shape[1])
        
        fold += 1
    
    # plt.plot(train_loss, label="Training loss")
    # plt.plot(val_loss, label="Validation Loss")
    # plt.title("Training and Validation Loss")
    # plt.xlabel("Epochs")
    # plt.ylabel("Loss")
    # plt.legend()
    # plt.savefig(logging_info["folder_path"] + "plots/" + logging_info["file_name"] + "_loss.png")
    # plt.close()

    # plt.plot(val_corr, label="Validation Correlation")
    # plt.title("Validation Correlation")
    # plt.xlabel("Epochs")
    # plt.ylabel("Correlation")
    # plt.savefig(logging_info["folder_path"] + "plots/" + logging_info["file_name"] + "_correlation.png")
