from train import *
from datasets import *
import sys


hyperparameters = {"conv_hidden_channels":64, "mlp_hidden_channels":256, "num_epochs":200, "learning_rate":0.0001, "batch_size":32,
                   "mlp_num_layers":2, "edge_lengths":[[1,2,3],[5,10,15,20,25]], "kernel_length":2,
                   "max_time":8*60*60, "conv_dropout": 0.2, "mlp_dropout": 0.2}

tiling_data = DNASeqDataset(src_file="~/Dropbox/Northwestern/DNA_Cyclizability/cycle5.txt", 
                            edge_lengths=hyperparameters["edge_lengths"], kernel_length=hyperparameters["kernel_length"])
random_data = DNASeqDataset(src_file="~/Dropbox/Northwestern/DNA_Cyclizability/cycle3.txt", 
                            edge_lengths=hyperparameters["edge_lengths"], kernel_length=hyperparameters["kernel_length"])

logging_info = {"print_every":1, "folder_path":"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/benchmarks/cnns/",
                "file_name": "edge_separate_1_2_3__5_10_15_20_25_with_mlps_3layer_kernellength2_64convhc_256mlphc_02convdropout_02mlpdropout_8hr"}
log_file_path = logging_info["folder_path"] + "log_files/" + logging_info["file_name"] + ".txt"

with open(log_file_path, "w") as f:
    sys.stdout = f

    print(f"Hyperparameters: \n \t {hyperparameters}")

    print(f"\n \n \n Begin Training: \n \n \n")

    model, train_loss, val_loss, val_corr = train(tiling_data, hyperparameters=hyperparameters, 
                                                  logging_info=logging_info, 
                                                  model_name="ConvSeparateEdgeLengthsWithMLPs", 
                                                  val_data=random_data, plot_0_val=False)
    
    plt.plot(train_loss, label="Training loss")
    plt.plot(val_loss, label="Validation Loss")
    plt.title("Training and Validation Loss")
    plt.xlabel("Epochs")
    plt.ylabel("Loss")
    plt.legend()
    plt.savefig(logging_info["folder_path"] + "plots/" + logging_info["file_name"] + "_loss.png")
    plt.close()

    plt.plot(val_corr, label="Validation Correlation")
    plt.title("Validation Correlation")
    plt.xlabel("Epochs")
    plt.ylabel("Correlation")
    plt.savefig(logging_info["folder_path"] + "plots/" + logging_info["file_name"] + "_correlation.png")
