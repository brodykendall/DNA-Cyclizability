from train import *
# from process_data import *
from datasets import *
import sys

# CustomUNet:

# tiling_data = DNASeqGraph("/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/benchmarks/gnns", 
#                           raw_graph_name="cycle5", processed_graph_name="cycle5_1_to_49_edge_together_include_indices", 
#                           edge_distances=list(range(1,50)), include_indices=True)
# random_data = DNASeqGraph("/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/benchmarks/gnns", 
#                           raw_graph_name="cycle3", processed_graph_name="cycle3_1_to_49_edge_together_include_indices", 
#                           edge_distances=list(range(1,50)), include_indices=True)
# 
# cluster_list = [
#     torch.Tensor(
#         [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])
# ]
# hyperparameters = {"hidden_channels":1024, "num_layers":[1], "num_epochs":200, "learning_rate":0.0001, "batch_size":32,
#                    "pool_ratios":[], "pool_types":[], "mlp_num_layers":2, "edge_lengths":list(range(1,50)), 
#                    "mp_layer":GAT, "jk":"lstm", "edge_weights":False, "cluster_list":cluster_list, "max_time":18*60*60}


# GraphSAGESeparateEdgeLengths:

# tiling_data = DNASeqGraphSeparateEdges("/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/benchmarks/gnns", 
#                                        raw_graph_name="cycle5", 
#                                        processed_graph_name="cycle5_1_to_49_edge_separate_include_indices", 
#                                        edge_distances=list(range(1,50)),
#                                        include_indices=True)
# random_data = DNASeqGraphSeparateEdges("/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/benchmarks/gnns", 
#                                        raw_graph_name="cycle3", 
#                                        processed_graph_name="cycle3_1_to_49_edge_separate_include_indices", 
#                                        edge_distances=list(range(1,50)),
#                                        include_indices=True)
# 
# cluster_list = [
#     torch.Tensor(
#         [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])
# ]
# hyperparameters = {"hidden_channels":128, "num_layers":[2], "num_epochs":200, "learning_rate":0.0001, "batch_size":32,
#                    "pool_ratios":[], "pool_types":[], "mlp_num_layers":2, "edge_lengths":[1,5,10,15], 
#                    "mp_layer":GCN, "jk":"lstm", "cluster_list":cluster_list, "max_time":18*60*60}


# GraphSAGESeparateEdgeLengthsWithMLPs:

train_data = DNASeqGraphSeparateEdges("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/benchmarks/gnns", 
                                       raw_graph_dir = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/",
                                       raw_graph_filetype="csv",
                                       raw_graph_name="yeast_chrV_ir_lstm_tiling_post_smoothed_matched", 
                                       variable_of_interest="model_ir_lstm_tiling_pred_chrV_post_smooth_5bp",
                                       processed_graph_name="post_smoothed_1_2_3_5_10_15_20_25_30_edge_separate_include_indices", 
                                       edge_distances=[1,2,3,5,10,15,20,25,30],
                                       include_indices=True)
# random_data = DNASeqGraphSeparateEdges("/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/benchmarks/gnns", 
#                                        raw_graph_name="cycle3", 
#                                        processed_graph_name="cycle3_1_2_3_5_10_15_20_25_30_edge_separate_include_indices", 
#                                        edge_distances=[1,2,3,5,10,15,20,25,30],
#                                        include_indices=True)

# cluster_list = [
#     torch.Tensor(
#         [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])
# ]
# hyperparameters = {"mp_hidden_channels":64, "mlp_hidden_channels":64, "num_layers":[1], "num_epochs":200, "learning_rate":0.0001, "batch_size":32,
#                    "pool_ratios":[], "pool_types":[], "mlp_num_layers":2, "edge_lengths":[[1,2,3],[5,10,15,20,25,30]], 
#                    "mp_layer":GCN, "jk":"lstm", "cluster_list":cluster_list, "max_time":4*60*60}

hyperparameters = {"mp_hidden_channels":48, "mlp_hidden_channels":128, "num_layers":[1], "num_epochs":200, "learning_rate":0.0001, "batch_size":32,
                   "pool_ratios":[], "pool_types":[], "mlp_num_layers":2, "edge_lengths":[[1,2,3],[5,10,15,20,25,30]], 
                   "mp_layer":GCN, "jk":"lstm", "num_nodes":50, "max_time":8*60*60, "final_dropout": 0.2}


logging_info = {"print_every":1, "folder_path":"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/benchmarks/gnns/",
                "file_name": "post_smoothed_edge_separate_1_2_3__5_10_15_20_25_30_with_mlps_2layer_GCN_1layer_48mphc_128mlphc_include_indices_no_final_pooling_02finaldropout_8hr"}
log_file_path = logging_info["folder_path"] + "log_files/" + logging_info["file_name"] + ".txt"

with open(log_file_path, "w") as f:
    sys.stdout = f

    print(f"Hyperparameters: \n \t {hyperparameters}")

    print(f"\n \n \n Begin Training: \n \n \n")

    # model, train_loss, val_loss, val_corr = train(train_data, hyperparameters=hyperparameters, 
    #                                               logging_info=logging_info, 
    #                                               model_name="GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling", 
    #                                               val_data=random_data, plot_0_val=False)

    model, train_loss, val_loss, val_corr = train(train_data, hyperparameters=hyperparameters, 
                                                  logging_info=logging_info, 
                                                  model_name="GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling", 
                                                  plot_0_val=False)
    
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
