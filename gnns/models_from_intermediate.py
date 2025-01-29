from imports import *

# GraphSAGESeparateEdgeLengthsWithMLPs:
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