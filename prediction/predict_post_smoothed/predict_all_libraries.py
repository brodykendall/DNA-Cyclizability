from models import *
from train import initialize_model

hyperparameters = {"nodes_per_layer":16,"num_epochs":200, "learning_rate":0.00001, "batch_size":32,
                   "num_layers":5, "max_time":8*60*60, "dropout": 0.0, "patience": 30}

logging_info = {"print_every":1, "folder_path":"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/predict_post_smoothed/",
                "file_name": f"5bp_DNAcycP_{hyperparameters['nodes_per_layer']}nodes_per_layer_{hyperparameters['num_layers']}layers_{int(hyperparameters['dropout']*10)}dropout_3input_features"}
# logging_info = {"print_every":1, "folder_path":"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/predict_post_smoothed/",
#                 "file_name": f"{hyperparameters['nodes_per_layer']}nodes_per_layer_{hyperparameters['num_layers']}layers_{hyperparameters['dropout']}dropout"}
# logging_info = {"print_every":1, "folder_path":"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/predict_post_smoothed/",
#                 "file_name": f"{hyperparameters['nodes_per_layer']}nodes_per_layer_{hyperparameters['num_layers']}layers_{hyperparameters['dropout']}dropout_inc_fourier"}

model = initialize_model("FeedForward", 1, hyperparameters, 3)
model.load_state_dict(torch.load(logging_info["folder_path"] + "models/fold0/" + logging_info["file_name"]))

# cycle1 = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/cycle1.txt",delimiter = ",")
# X1 = torch.Tensor(np.array(cycle1[["n=26", "n=29", "n=31"]]))
# with torch.no_grad():
#     X1_pred = model(X1)
# cycle1["predicted_smooth_C0"] = X1_pred
# cycle1.to_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/cycle1_predicted_smooth_C0.csv", index=False)

cycle3 = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/cycle3.txt",delimiter = ",")
X3 = torch.Tensor(np.array(cycle3[["n=26", "n=29", "n=31"]]))
with torch.no_grad():
    X3_pred = model(X3)
cycle3["predicted_smooth_C0"] = X3_pred
cycle3.to_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/cycle3_predicted_smooth_C0.csv", index=False)

cycle5 = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/cycle5.txt",delimiter = ",")
X5= torch.Tensor(np.array(cycle5[["n=26", "n=29", "n=31"]]))
with torch.no_grad():
    X5_pred = model(X5)
cycle5["predicted_smooth_C0"] = X5_pred
cycle5.to_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/cycle5_predicted_smooth_C0.csv", index=False)

cycle6 = pd.read_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/cycle6.txt",delimiter = ",")
X6 = torch.Tensor(np.array(cycle6[["n=26", "n=29", "n=31"]]))
with torch.no_grad():
    X6_pred = model(X6)
cycle6["predicted_smooth_C0"] = X6_pred
cycle6.to_csv("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/cycle6_predicted_smooth_C0.csv", index=False)
