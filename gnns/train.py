from imports import *
from helper_functions import *
from models import *
import time
import networkx as nx

def initialize_model(model_name:str, num_node_features, output_features, hyperparameters):
    if model_name == "CustomUNet":
        net = CustomUNet(
            num_node_features=num_node_features, hidden_channels=hyperparameters["hidden_channels"], 
            num_layers=hyperparameters["num_layers"], output_features=output_features, 
            pool_ratios=hyperparameters["pool_ratios"], pool_types=hyperparameters["pool_types"], 
            mlp_num_layers=hyperparameters["mlp_num_layers"], jk=hyperparameters["jk"], 
            cluster_list = hyperparameters["cluster_list"], mp_layer=hyperparameters["mp_layer"], 
            edge_weights=hyperparameters["edge_weights"])
    elif model_name == "GraphSAGESeparateEdgeLengths":
        cluster_list=hyperparameters["cluster_list"] if "cluster_list" in hyperparameters.keys() else []
        net = GraphSAGESeparateEdgeLengths(
            num_node_features=num_node_features, hidden_channels=hyperparameters["hidden_channels"], 
            num_layers=hyperparameters["num_layers"], output_features=output_features, 
            pool_ratios=hyperparameters["pool_ratios"], pool_types=hyperparameters["pool_types"], 
            mlp_num_layers=hyperparameters["mlp_num_layers"], edge_lengths=hyperparameters["edge_lengths"], 
            mp_layer=hyperparameters["mp_layer"], jk=hyperparameters["jk"], cluster_list=cluster_list)
    elif model_name == "HierarchicalPoolingSeparateEdgeLengths":
        net = HierarchicalPoolingSeparateEdgeLengths(
            num_node_features=num_node_features, hidden_channels=hyperparameters["hidden_channels"], 
            num_layers=hyperparameters["num_layers"], output_features=output_features, 
            pool_ratios=hyperparameters["pool_ratios"], pool_types=hyperparameters["pool_types"], 
            mlp_num_layers=hyperparameters["mlp_num_layers"], edge_lengths=hyperparameters["edge_lengths"], 
            jk=hyperparameters["jk"], cluster_list=hyperparameters["cluster_list"], 
            mp_layer=hyperparameters["mp_layer"])
    elif model_name == "GraphSAGESeparateEdgeLengthsWithMLPs":
        cluster_list=hyperparameters["cluster_list"] if "cluster_list" in hyperparameters.keys() else []
        net = GraphSAGESeparateEdgeLengthsWithMLPs(
            num_node_features=num_node_features, mp_hidden_channels=hyperparameters["mp_hidden_channels"], 
            mlp_hidden_channels=hyperparameters["mlp_hidden_channels"], num_layers=hyperparameters["num_layers"], 
            output_features=output_features, pool_ratios=hyperparameters["pool_ratios"], pool_types=hyperparameters["pool_types"], 
            mlp_num_layers=hyperparameters["mlp_num_layers"], edge_lengths=hyperparameters["edge_lengths"], 
            mp_layer=hyperparameters["mp_layer"], jk=hyperparameters["jk"], cluster_list=cluster_list)
    elif model_name == "GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling":
        net = GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling(
            num_node_features=num_node_features, mp_hidden_channels=hyperparameters["mp_hidden_channels"], 
            mlp_hidden_channels=hyperparameters["mlp_hidden_channels"], num_layers=hyperparameters["num_layers"], 
            output_features=output_features, pool_ratios=hyperparameters["pool_ratios"], pool_types=hyperparameters["pool_types"], 
            mlp_num_layers=hyperparameters["mlp_num_layers"], edge_lengths=hyperparameters["edge_lengths"], 
            mp_layer=hyperparameters["mp_layer"], jk=hyperparameters["jk"], num_nodes=hyperparameters["num_nodes"],
            final_dropout=hyperparameters["final_dropout"])
    return net

def train(train_data, hyperparameters, logging_info, model_name:str, val_data=None, test_data=None, plot_0_val=False):
    num_epochs = hyperparameters["num_epochs"] if "num_epochs" in hyperparameters.keys() else 20
    learning_rate = hyperparameters["learning_rate"] if "learning_rate" in hyperparameters.keys() else 0.01
    batch_size = hyperparameters["batch_size"] if "batch_size" in hyperparameters.keys() else 20
    lr_reduce_threshold = hyperparameters["lr_reduce_threshold"] if "lr_reduce_threshold" in hyperparameters.keys() else 0
    max_time = hyperparameters["max_time"] if "max_time" in hyperparameters.keys() else 4*60*60

    print_every = logging_info["print_every"] if "print_every" in logging_info.keys() else 5
    folder_path = logging_info["folder_path"]
    file_name = logging_info["file_name"]
    patience = logging_info["patience"] if "patience" in logging_info.keys() else 10

    model_save_path = folder_path + "/models/" + file_name

    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    if val_data:
        val_loader = DataLoader(val_data, batch_size=1)
    
    net = initialize_model(model_name=model_name, num_node_features=train_data[0].x.shape[1], 
                           output_features=train_data[0].y.shape[0], hyperparameters=hyperparameters)
    
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(net.parameters(), lr=learning_rate)
    scheduler = ReduceLROnPlateau(optimizer, 'min', patience=patience//2, 
                                  threshold=lr_reduce_threshold, verbose=True)

    training_loss = []
    validation_loss = []
    validation_corr = []

    best_loss = 1e9
    no_improvement_count = 0

    start_time = time.time()
    for epoch in range(num_epochs):
        
        num_batches=0
        train_loss = 0

        net.train(True)

        for _, train_batch in enumerate(train_loader):
            optimizer.zero_grad()
            pred = net(train_batch)
            batch_loss = criterion(pred, train_batch.y)
            batch_loss.backward()
            optimizer.step()
            
            
            num_batches += 1
            train_loss += batch_loss.item()
            
        training_loss.append(train_loss/num_batches)

        net.eval()

        with torch.no_grad():
            if epoch % print_every == 0:
                print(f"End of epoch: {epoch}. Total run time: {round(time.time() - start_time)} seconds")
                print(f"\t Train loss: {training_loss[-1]}", flush=True)
                if val_data:
                    num_val_batch = 0
                    val_loss = 0
                    val_predictions = []
                    val_true = []
                    for val_idx, val_batch in enumerate(val_loader):
                        val_pred = net(val_batch, plot=True) if plot_0_val and val_idx == 0 else net(val_batch)
                        val_batch_loss = criterion(val_pred, val_batch.y)

                        if np.random.rand() < 5/len(val_data):
                            print(f"\t Randomly selected predicted value: {val_pred}, true value: {val_batch.y}")

                        num_val_batch += 1
                        val_loss += val_batch_loss.item()
                        val_predictions.append(val_pred.item())
                        val_true.append(val_batch.y.item())
                    
                    validation_loss.append(val_loss/num_val_batch)
                    validation_corr.append(np.corrcoef((val_predictions, val_true))[0,1])

                    print(f"\t Validation loss: {validation_loss[-1]}")
                    print(f"\t Validation Correlation: {validation_corr[-1]}", flush=True)

            cur_loss = validation_loss[-1] if val_data else training_loss[-1]
            improved, best_loss, no_improvement_count = model_improved(cur_loss, best_loss, no_improvement_count, 
                                                                       net, model_save_path)
            if no_improvement_count >= patience:
                print(f"Model has not improved over the last {patience} epochs. Stopping early.")
                break
            if time.time() - start_time > max_time:
                print(f"Total time of {round(time.time() - start_time)} seconds has exceeded the maximum time: {max_time}. Stopping early.")
                break
            scheduler.step(cur_loss)

    net.load_state_dict(torch.load(model_save_path))
    return net, training_loss, validation_loss, validation_corr


        

