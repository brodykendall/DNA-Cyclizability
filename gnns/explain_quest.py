from imports import *
# from process_data import *
from train import initialize_model
from datasets import *
import matplotlib.pyplot as plt
import time
from torch_geometric.explain import Explainer
from torch_geometric.explain import ModelConfig
from torch_geometric.explain import ExplainerConfig
from torch_geometric.explain import DummyExplainer
from explain import *

import types
import networkx as nx

tiling_data = DNASeqGraphSeparateEdges("~/Documents/DNA_Cyclizability/benchmarks/gnns", 
                                       raw_graph_name="cycle5", 
                                       processed_graph_name="cycle5_1_2_3_5_10_15_20_25_30_edge_separate_include_indices", 
                                       edge_distances=[1,2,3,5,10,15,20,25,30],
                                       include_indices=True)
random_data = DNASeqGraphSeparateEdges("~/Documents/DNA_Cyclizability/benchmarks/gnns", 
                                       raw_graph_name="cycle3", 
                                       processed_graph_name="cycle3_1_2_3_5_10_15_20_25_30_edge_separate_include_indices", 
                                       edge_distances=[1,2,3,5,10,15,20,25,30],
                                       include_indices=True)
yeast_data = DNASeqGraphSeparateEdges("~/Documents/DNA_Cyclizability/benchmarks/gnns", 
                                      raw_graph_name="cycle1", 
                                      processed_graph_name="cycle1_1_2_3_5_10_15_20_25_30_edge_separate_include_indices", 
                                      edge_distances=[1,2,3,5,10,15,20,25,30],
                                      include_indices=True)
chrv_data = DNASeqGraphSeparateEdges("~/Documents/DNA_Cyclizability/benchmarks/gnns", 
                                     raw_graph_name="cycle6", 
                                     processed_graph_name="cycle6_1_2_3_5_10_15_20_25_30_edge_separate_include_indices", 
                                     edge_distances=[1,2,3,5,10,15,20,25,30],
                                     include_indices=True)


hyperparameters = {"mp_hidden_channels":48, "mlp_hidden_channels":128, "num_layers":[1], "num_epochs":200, "learning_rate":0.0001, "batch_size":32,
                   "pool_ratios":[], "pool_types":[], "mlp_num_layers":2, "edge_lengths":[[1,2,3],[5,10,15,20,25,30]], 
                   "mp_layer":GCN, "jk":"lstm", "num_nodes":50, "max_time":8*60*60, "final_dropout": 0.2}

data_str = "random"

if data_str == "tiling":
    data = tiling_data
elif data_str == "random":
    data = random_data
elif data_str == "yeast":
    data = yeast_data
elif data_str == "chrv":
    data = chrv_data

num_obs=len(data)

logging_info = {"print_every":1, "folder_path":"/home/cbk6686/Documents/DNA_Cyclizability/benchmarks/gnns/",
                "file_name": "edge_separate_1_2_3__5_10_15_20_25_30_with_mlps_2layer_GCN_1layer_48mphc_128mlphc_include_indices_no_final_pooling_02finaldropout_8hr"}

model = initialize_model(model_name="GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling", num_node_features=data.x.shape[1],
                         output_features=data.y.shape[1], hyperparameters=hyperparameters)
model.load_state_dict(torch.load(logging_info["folder_path"] + "models/" + logging_info["file_name"]))
model.eval()



_, which_graph = torch.topk(data.y.flatten(), k=num_obs)
which_graph = np.array(which_graph)
explanation_data = DataLoader(data[which_graph], batch_size=1)


def forward_until_mlps(model, data):
    x, batch = data.x, data.batch
    network_outputs = []

    network_idx = 0

    for outer_idx, edge_length_list in enumerate(model.edge_lengths):
        inner_loop_network_outputs = []
        for edge_length in edge_length_list:

            # edge_length_attr_name = f"edge_index_{edge_length}_locations"
            edge_index_locations = abs(data.edge_index[0,:] - data.edge_index[1,:])==edge_length
            edge_index_locations += (data.edge_index[0,:] - data.edge_index[1,:]) == 0
            # edge_index_locations = getattr(data, edge_length_attr_name, None)
            edge_index = data.edge_index[:,edge_index_locations]
            if edge_index is None:
                raise ValueError(f"Edge length of {edge_length} not found in data")
            
            current_network = model.edge_length_networks[network_idx]

            x_cur = current_network.convs[0](x=x, edge_index=edge_index, edge_weight=None)
            batch_cur = batch

            x_pooling = []
            x_pooling.append(x_cur)

            pooling_indices = []
            pooling_edge_indices = []

            for i, pooler in enumerate(current_network.poolers):
                if current_network.pool_types[i] == "ASA":
                    x_cur, edge_index, edge_weight, batch_cur, index = pooler(
                        x=x_cur, edge_index=edge_index, batch=batch_cur)
                else:
                    x_cur, edge_index, edge_weight, batch_cur, index, score = pooler(
                        x=x_cur, edge_index=edge_index, batch=batch_cur)
                
                pooling_indices.append(index)
                pooling_edge_indices.append(edge_index)
                
                x_cur = current_network.convs[i+1](x=x_cur, edge_index=edge_index)
                if i < len(current_network.poolers) - 1:
                    x_pooling.append(x_cur)
            
            for i, unpooler in enumerate(current_network.unpoolers):
                unpooling_indices = pooling_indices.pop()
                unpooling_edge_index = pooling_edge_indices.pop()
                x_unpooling = torch.zeros(x_pooling[-1].shape)
                x_unpooling[unpooling_indices] = x_cur
                x_unpooling = unpooler(x=x_unpooling, edge_index=unpooling_edge_index)

                x_unpooling = torch.cat((x_unpooling, x_pooling.pop()), dim=1)
                
                x_cur = current_network.unpooling_mlps[i](x_unpooling)

            inner_loop_network_outputs.append(current_network.final_mlp(x_cur))

            network_idx += 1
        
        network_outputs.append(torch.cat(inner_loop_network_outputs, dim=1))

    return torch.cat(network_outputs, dim=1)

def edge_length_final_mlps(x, model, hidden_channels):
    outer_loop_network_outputs = []
    prev_idx = 0
    for outer_idx, edge_length_list in enumerate(model.edge_lengths):
        x_cur = x[:,prev_idx:(prev_idx+hidden_channels*len(edge_length_list))]
        outer_loop_network_outputs.append(model.edge_lengths_final_mlps[outer_idx](x_cur))
        prev_idx += hidden_channels*len(edge_length_list)
    return torch.cat(outer_loop_network_outputs, dim=1)

def final_mlps_no_pooling(x, model, batch):
    x = model.final_mlp_1(x)
    x = x.reshape(torch.max(batch)+1, -1)
    x = model.final_mlp_2(x)
    return x

def explain_message(self, inputs: torch.Tensor, size_i: int) -> torch.Tensor:
        # NOTE Replace this method in custom explainers per message-passing
        # layer to customize how messages shall be explained, e.g., via:
        # conv.explain_message = explain_message.__get__(conv, MessagePassing)
        # see stackoverflow.com: 394770/override-a-method-at-instance-level

    edge_mask = self._edge_mask

    if edge_mask is None:
        raise ValueError(f"Could not find a pre-defined 'edge_mask' as "
                         f"part of {self.__class__.__name__}.")

    if self._apply_sigmoid:
        edge_mask = edge_mask.sigmoid()

    # batch_size = num_obs//10
    batch_size = 1
    if inputs.size(self.node_dim) == 148*batch_size:
        # Edge length 1 (98)
        edge_mask = edge_mask[0:98*batch_size]
    if inputs.size(self.node_dim) == 146*batch_size:
        # Edge length 2 (96)
        edge_mask = edge_mask[98*batch_size:194*batch_size]
    if inputs.size(self.node_dim) == 144*batch_size:
        # Edge length 3 (94)
        edge_mask = edge_mask[194*batch_size:288*batch_size]
    elif inputs.size(self.node_dim) == 140*batch_size:
        # Edge length 5 (90)
        edge_mask = edge_mask[288*batch_size:378*batch_size]
    elif inputs.size(self.node_dim) == 130*batch_size:
        # Edge length 10 (80)
        edge_mask = edge_mask[378*batch_size:458*batch_size]
    elif inputs.size(self.node_dim) == 120*batch_size:
        # Edge length 15 (70)
        edge_mask = edge_mask[458*batch_size:528*batch_size]
    if inputs.size(self.node_dim) == 110*batch_size:
        # Edge length 20 (60)
        edge_mask = edge_mask[528*batch_size:588*batch_size]
    if inputs.size(self.node_dim) == 100*batch_size:
        # Edge length 25 (50)
        edge_mask = edge_mask[588*batch_size:638*batch_size]
    if inputs.size(self.node_dim) == 90*batch_size:
        # Edge length 30 (40)
        edge_mask = edge_mask[638*batch_size:678*batch_size]

     # Some ops add self-loops to `edge_index`. We need to do the same for
     # `edge_mask` (but do not train these entries).
    if inputs.size(self.node_dim) != edge_mask.size(0):
        # edge_mask = edge_mask[self._loop_mask]
        loop = edge_mask.new_ones(size_i)
        edge_mask = torch.cat([edge_mask, loop], dim=0)
    if inputs.size(self.node_dim) != edge_mask.size(0):
        print(f"\t Inputs: {inputs.shape} \n \t Edge mask: {edge_mask.shape}")
        assert inputs.size(self.node_dim) == edge_mask.size(0)

    size = [1] * inputs.dim()
    size[self.node_dim] = -1
    return inputs * edge_mask.view(size)

for i, module in enumerate(model.modules()):
    if hasattr(module, "explain_message"):
        # funcType = type(module.explain_message)
        # module.explain_message = funcType(explain_message2, module)
        module.explain_message = types.MethodType(explain_message, module)

model_config = ModelConfig(mode="regression", task_level="graph", return_type="raw")
explainer_config = ExplainerConfig(explanation_type="model", edge_mask_type="object")

model_2 = tg_nn.Sequential('x, batch, orig_model, hidden_channels', [
    (edge_length_final_mlps, 'x, orig_model, hidden_channels -> x'), 
    (final_mlps_no_pooling, 'x, orig_model, batch -> x')])

test_explainer = CustomExplainer(lr = 0.1)
test_explainer.connect(explainer_config=explainer_config, model_config=model_config)
test_explainer.epochs=100

test_explanations=[]
for i, explanation_batch in enumerate(explanation_data):
    test_explanations.append(test_explainer(model=model, data=explanation_batch, target=explanation_batch.y))
    if (i+1) % 1000 == 0:
        torch.save(test_explanations, logging_info["folder_path"] + "explanations/" + logging_info["file_name"] + "_" + data_str + "_full")
        print(f"Completed batch {i}")

torch.save(test_explanations, logging_info["folder_path"] + "explanations/" + logging_info["file_name"] + "_" + data_str + "_full")