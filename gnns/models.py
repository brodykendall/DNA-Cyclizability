from imports import *
import matplotlib.pyplot as plt
import networkx as nx


class CustomUNet(Module):
    def __init__(self, num_node_features, mp_hidden_channels, mlp_hidden_channels, num_layers, output_features, pool_ratios, pool_types, 
                 mlp_num_layers, jk:str, cluster_list:list[torch.Tensor]=[], mp_layer=GraphSAGE, edge_weights:bool=False):
        super(CustomUNet, self).__init__()

        self.convs = ModuleList()
        self.poolers = ModuleList()
        self.unpoolers = ModuleList()
        self.unpooling_mlps = ModuleList()
        self.pool_types = pool_types
        self.edge_weights=edge_weights
        self.cluster_list = cluster_list

        self.convs.append(mp_layer(in_channels=num_node_features, hidden_channels=mp_hidden_channels, 
            num_layers=num_layers[0], output_features=mp_hidden_channels, jk=jk))
        for i, pool_ratio in enumerate(pool_ratios):
            if self.pool_types[i] == "ASA":
                self.poolers.append(ASAPooling(in_channels=mp_hidden_channels, ratio=pool_ratio))
            else:
                self.poolers.append(SAGPooling(in_channels=mp_hidden_channels, ratio=pool_ratio))
            self.convs.append(mp_layer(in_channels=mp_hidden_channels, hidden_channels=mp_hidden_channels, 
                                        num_layers=num_layers[i+1], output_features=mp_hidden_channels, jk=jk))
            self.unpoolers.append(mp_layer(in_channels=mp_hidden_channels, hidden_channels=mlp_hidden_channels,
                                            num_layers=num_layers[i+1], output_features=mlp_hidden_channels, jk=jk))
            self.unpooling_mlps.append(MLP(in_channels=mlp_hidden_channels*2, hidden_channels=mlp_hidden_channels, 
                                           num_layers=mlp_num_layers, out_channels=mlp_hidden_channels))
        # self.final_mlp_1 = MLP(in_channels=hidden_channels, hidden_channels=hidden_channels,
        #     num_layers=mlp_num_layers, out_channels=hidden_channels)
        self.final_mlp = MLP(in_channels=mp_hidden_channels, hidden_channels=mlp_hidden_channels,
            num_layers=mlp_num_layers, out_channels=output_features)
    
    def forward(self, data, plot=False):
    # def forward(self, data):
        plot=False
        x, edge_index, batch = data.x, data.edge_index, data.batch
        if self.edge_weights:
            edge_weight = data.edge_weight
            x = self.convs[0](x=x, edge_index=edge_index, edge_weight=edge_weight)
        else:
            x = self.convs[0](x=x, edge_index=edge_index)

        x_pooling = []
        x_pooling.append(x)

        pooling_indices = []
        pooling_edge_indices = []
        if self.edge_weights:
            pooling_edge_weights = []

        if plot:
            plt.figure(figsize=(14,12))
            g = torch_geometric.utils.to_networkx(data, to_undirected=True, remove_self_loops=True)
            plotting_labels = {node: str(node) for node in g.nodes()}
            plotting_pos = nx.circular_layout(g)
            nx.draw(g, labels=plotting_labels, pos=plotting_pos)
            plt.title("Original")
            plt.show()

        for i, pooler in enumerate(self.poolers):
            if self.pool_types[i] == "ASA":
                if self.edge_weights:
                    x, edge_index, edge_weight, batch, index = pooler(
                        x=x, edge_index=edge_index, batch=batch, edge_weight=edge_weight)
                else:
                    x, edge_index, edge_weight, batch, index = pooler(
                        x=x, edge_index=edge_index, batch=batch)
            else:
                if self.edge_weights:
                    x, edge_index, edge_weight, batch, index, score = pooler(
                        x=x, edge_index=edge_index, batch=batch, edge_weight=edge_weight)
                else:
                    x, edge_index, edge_weight, batch, index, score = pooler(
                        x=x, edge_index=edge_index, batch=batch)
            
            pooling_indices.append(index)
            pooling_edge_indices.append(edge_index)
            if self.edge_weights:
                pooling_edge_weights.append(edge_weight)
            
            if plot:
                plt.figure(figsize=(14,12))
                plot_data = Data(x=x, edge_index=edge_index)
                g = torch_geometric.utils.to_networkx(plot_data, to_undirected=True, remove_self_loops=True)
                plotting_labels = {node: plotting_labels[index[node].item()] for node in g.nodes()}
                plotting_pos = {node: plotting_pos[index[node].item()] for node in g.nodes()}
                nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                plt.title(f"Pooled level {i}")
                plt.show()
            
            if self.edge_weights:
                x = self.convs[i+1](x=x, edge_index=edge_index, edge_weight=edge_weight)
            else:
                x = self.convs[i+1](x=x, edge_index=edge_index)
            
            if i < len(self.poolers) - 1:
                x_pooling.append(x)
        
        for i, unpooler in enumerate(self.unpoolers):
            unpooling_indices = pooling_indices.pop()
            unpooling_edge_index = pooling_edge_indices.pop()
            x_unpooling = torch.zeros(x_pooling[-1].shape)
            x_unpooling[unpooling_indices] = x
            if self.edge_weights:
                unpooling_edge_weight = pooling_edge_weights.pop()
                x_unpooling = unpooler(x=x_unpooling, edge_index=unpooling_edge_index, edge_weight=unpooling_edge_weight)
            else:
                x_unpooling = unpooler(x=x_unpooling, edge_index=unpooling_edge_index)

            x_unpooling = torch.cat((x_unpooling, x_pooling.pop()), dim=1)
            
            x = self.unpooling_mlps[i](x_unpooling)

        batch = data.batch
        # x = self.final_mlp_1(x)
        for cluster in self.cluster_list:
            cur_cluster = cluster.repeat(len(batch)//len(cluster)) + (torch.arange(len(batch))//len(cluster))*(torch.max(cluster)+1)
            x, batch = max_pool_x(cur_cluster, x, batch)
        x = global_mean_pool(x, batch)
        # x = global_mean_pool(x, data.batch)
        # x = global_add_pool(x, data.batch)

        x = self.final_mlp(x)

        return x
    
class GraphSAGESeparateEdgeLengths(Module):
    def __init__(self, num_node_features, hidden_channels, num_layers, output_features, pool_ratios, pool_types, 
                 mlp_num_layers, edge_lengths, jk:str, mp_layer=GraphSAGE, cluster_list:list=[]):
        super(GraphSAGESeparateEdgeLengths, self).__init__()
        self.edge_lengths = edge_lengths
        self.edge_length_networks = ModuleList()
        self.cluster_list = cluster_list

        for _ in edge_lengths:
            self.edge_length_networks.append(CustomUNet(num_node_features=num_node_features, 
                                                        hidden_channels=hidden_channels,
                                                        num_layers=num_layers,
                                                        output_features=hidden_channels,
                                                        pool_ratios=pool_ratios,
                                                        pool_types=pool_types,
                                                        mlp_num_layers=mlp_num_layers,
                                                        jk=jk,
                                                        mp_layer=mp_layer))
            
        self.final_mlp_1 = MLP(in_channels=hidden_channels*len(edge_lengths), hidden_channels=hidden_channels,
                               num_layers=mlp_num_layers, out_channels=hidden_channels)
        self.final_mlp_2 = MLP(in_channels=hidden_channels, hidden_channels=hidden_channels,
                               num_layers=mlp_num_layers, out_channels=output_features)
    
    def forward(self, data, plot=False):
    # def forward(self, data):
        # plot=False
        x, batch = data.x, data.batch
        network_outputs = []

        for network_idx, edge_length in enumerate(self.edge_lengths):

            # edge_length_attr_name = f"edge_index_{edge_length}_locations"
            edge_index_locations = abs(data.edge_index[0,:] - data.edge_index[1,:])==edge_length
            edge_index_locations += (data.edge_index[0,:] - data.edge_index[1,:]) == 0
            # edge_index_locations = getattr(data, edge_length_attr_name, None)
            edge_index = data.edge_index[:,edge_index_locations]
            if edge_index is None:
                raise ValueError(f"Edge length of {edge_length} not found in data")
            
            current_network = self.edge_length_networks[network_idx]

            x_cur = current_network.convs[0](x=x, edge_index=edge_index, edge_weight=None)
            batch_cur = batch

            x_pooling = []
            x_pooling.append(x_cur)

            pooling_indices = []
            pooling_edge_indices = []

            if plot:
                if network_idx%5 == 0:
                    cur_data = Data(x=x_cur, edge_index=edge_index)
                    plt.figure(figsize=(14,12))
                    g = torch_geometric.utils.to_networkx(cur_data, to_undirected=True, remove_self_loops=True)
                    plotting_labels = {node: str(node) for node in g.nodes()}
                    plotting_pos = nx.circular_layout(g)
                    nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                    plt.title(f"Original, edge length: {edge_length}")
                    plt.show()

            for i, pooler in enumerate(current_network.poolers):
                if current_network.pool_types[i] == "ASA":
                    x_cur, edge_index, edge_weight, batch_cur, index = pooler(
                        x=x_cur, edge_index=edge_index, batch=batch_cur)
                else:
                    x_cur, edge_index, edge_weight, batch_cur, index, score = pooler(
                        x=x_cur, edge_index=edge_index, batch=batch_cur)
                
                pooling_indices.append(index)
                pooling_edge_indices.append(edge_index)
                
                if plot:
                    if network_idx%5 == 0:
                        plt.figure(figsize=(14,12))
                        plot_data = Data(x=x_cur, edge_index=edge_index)
                        g = torch_geometric.utils.to_networkx(plot_data, to_undirected=True, remove_self_loops=True)
                        plotting_labels = {node: plotting_labels[index[node].item()] for node in g.nodes()}
                        plotting_pos = {node: plotting_pos[index[node].item()] for node in g.nodes()}
                        nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                        plt.title(f"Pooled level {i}, edge length: {edge_length}")
                        plt.show()
                
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

            network_outputs.append(current_network.final_mlp(x_cur))

        x = torch.cat(network_outputs, dim=1)
        x = self.final_mlp_1(x)
        batch = data.batch
        for cluster in self.cluster_list:
            cur_cluster = cluster.repeat(len(batch)//len(cluster)) + (torch.arange(len(batch))//len(cluster))*(torch.max(cluster)+1)
            x, batch = max_pool_x(cur_cluster, x, batch)
        x = global_mean_pool(x, batch)
        x = self.final_mlp_2(x)
        return x


class HierarchicalPoolingSeparateEdgeLengths(Module):
    def __init__(self, num_node_features, hidden_channels, num_layers, output_features, pool_ratios, pool_types, 
                 mlp_num_layers, edge_lengths, jk:str, cluster_list:list[torch.Tensor], mp_layer=GraphSAGE):
        super(HierarchicalPoolingSeparateEdgeLengths, self).__init__()
        self.edge_lengths = edge_lengths
        self.cluster_list = cluster_list
        # Note: If pooling_steps == 0, this should the same as GraphSAGESeparateEdgeLengths
        self.network_lists_by_pool_level = ModuleList()
        self.final_mlps_by_pool_level = ModuleList()

        # Top level Custom UNet:
        self.network_lists_by_pool_level.append([])
        for _ in edge_lengths:
            self.network_lists_by_pool_level[0].append(CustomUNet(
                num_node_features=num_node_features, hidden_channels=hidden_channels, 
                num_layers=num_layers, output_features=hidden_channels, 
                pool_ratios=pool_ratios, pool_types=pool_types, mlp_num_layers=mlp_num_layers, 
                jk=jk, mp_layer=mp_layer))
        self.final_mlps_by_pool_level.append(
            MLP(in_channels=hidden_channels*len(edge_lengths), hidden_channels=hidden_channels,
                num_layers=mlp_num_layers, out_channels=hidden_channels)
        )

        # Intermediate level Custom UNets:
        for i, _ in enumerate(cluster_list):
            self.network_lists_by_pool_level.append([])
            for _ in edge_lengths:
                self.network_lists_by_pool_level[i+1].append(CustomUNet(
                    num_node_features=hidden_channels, hidden_channels=hidden_channels, 
                    num_layers=num_layers, output_features=hidden_channels, 
                    pool_ratios=pool_ratios, pool_types=pool_types, mlp_num_layers=mlp_num_layers,
                    jk=jk, mp_layer=mp_layer))
            self.final_mlps_by_pool_level.append(
                MLP(in_channels=hidden_channels*len(edge_lengths), hidden_channels=hidden_channels,
                    num_layers=mlp_num_layers, out_channels=hidden_channels)
            )
        
        self.final_mlp = MLP(in_channels=hidden_channels, hidden_channels=hidden_channels,
                             num_layers=mlp_num_layers, out_channels=output_features)
    
    def forward(self, data, plot=False):
    # def forward(self, data):
        # plot=False
        for pool_level, network_list in enumerate(self.network_lists_by_pool_level):
            network_outputs = []
            if pool_level == 0:
                cur_level_data = data
            if pool_level != 0:
                cur_level_data = next_level_data
                cur_cluster = self.cluster_list[pool_level-1]
                cur_cluster = cur_cluster.repeat(len(cur_level_data.x)//len(cur_cluster)) + (torch.arange(len(cur_level_data.x))//len(cur_cluster))*(torch.max(cur_cluster)+1)
            
            next_level_data = Data()
            x_cur_level, batch_cur_level = cur_level_data.x, cur_level_data.batch
            
            for network_idx, edge_length in enumerate(self.edge_lengths):
                edge_length_attr_name = f"edge_index_{edge_length}"
                edge_index = getattr(cur_level_data, edge_length_attr_name, None)
                if edge_index is None:
                    raise ValueError(f"Edge length of {edge_length} not found in data")
                
                if pool_level == 0:
                    x_cur = x_cur_level
                    batch_cur = batch_cur_level

                else:
                    cur_data = Data(x=x_cur_level, edge_index=edge_index, batch=batch_cur_level)
                    cur_data = max_pool(cluster=cur_cluster, data=cur_data)
                    x_cur, edge_index, batch_cur = cur_data.x, cur_data.edge_index, cur_data.batch

                setattr(next_level_data, edge_length_attr_name, edge_index)
                
                current_network = network_list[network_idx]

                x_cur = current_network.convs[0](x=x_cur, edge_index=edge_index)

                x_pooling = []
                x_pooling.append(x_cur)

                pooling_indices = []
                pooling_edge_indices = []

                if plot:
                    if network_idx%5 == 0:
                        print("This plotting needs to be fixed to show the proper clusters")
                        cur_data = Data(x=x_cur, edge_index=edge_index)
                        plt.figure(figsize=(14,12))
                        g = torch_geometric.utils.to_networkx(cur_data, to_undirected=True, remove_self_loops=True)
                        plotting_labels = {node: str(node) for node in g.nodes()}
                        plotting_pos = nx.circular_layout(g)
                        nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                        plt.title(f"Original, edge length: {edge_length}")
                        plt.show()

                for i, pooler in enumerate(current_network.poolers):
                    if current_network.pool_types[i] == "ASA":
                        x_cur, edge_index, edge_weight, batch_cur, index = pooler(
                            x=x_cur, edge_index=edge_index, batch=batch_cur)
                    else:
                        x_cur, edge_index, edge_weight, batch_cur, index, score = pooler(
                            x=x_cur, edge_index=edge_index, batch=batch_cur)
                    
                    pooling_indices.append(index)
                    pooling_edge_indices.append(edge_index)
                    
                    if plot:
                        if network_idx%5 == 0:
                            plt.figure(figsize=(14,12))
                            plot_data = Data(x=x_cur, edge_index=edge_index)
                            g = torch_geometric.utils.to_networkx(plot_data, to_undirected=True, remove_self_loops=True)
                            plotting_labels = {node: plotting_labels[index[node].item()] for node in g.nodes()}
                            plotting_pos = {node: plotting_pos[index[node].item()] for node in g.nodes()}
                            nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                            plt.title(f"Pooled level {i}, edge length: {edge_length}")
                            plt.show()
                    
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

                network_outputs.append(current_network.final_mlp(x_cur))

            x = torch.cat(network_outputs, dim=1)
            x = self.final_mlps_by_pool_level[pool_level](x)
            next_level_data.x = x
            next_level_data.batch = batch_cur

        # x = torch.cat(network_outputs, dim=1)
        # x = self.final_mlp_1(x)
        # batch = batch_cur
        x = global_mean_pool(next_level_data.x, next_level_data.batch)
        x = self.final_mlp(x)
        return x

class GraphSAGESeparateEdgeLengthsWithMLPs(Module):
    def __init__(self, num_node_features, mp_hidden_channels, mlp_hidden_channels, num_layers, output_features, pool_ratios, pool_types, 
                 mlp_num_layers, edge_lengths, jk:str, mp_layer=GraphSAGE, cluster_list:list=[]):
        super(GraphSAGESeparateEdgeLengthsWithMLPs, self).__init__()
        self.edge_lengths = edge_lengths
        self.edge_length_networks = ModuleList()
        self.edge_lengths_final_mlps = ModuleList()
        self.cluster_list = cluster_list

        for edge_length_list in edge_lengths:
            for _ in edge_length_list:
                self.edge_length_networks.append(CustomUNet(num_node_features=num_node_features, 
                                                            mp_hidden_channels=mp_hidden_channels,
                                                            mlp_hidden_channels=mlp_hidden_channels,
                                                            num_layers=num_layers,
                                                            output_features=mlp_hidden_channels,
                                                            pool_ratios=pool_ratios,
                                                            pool_types=pool_types,
                                                            mlp_num_layers=mlp_num_layers,
                                                            jk=jk,
                                                            mp_layer=mp_layer))
            self.edge_lengths_final_mlps.append(MLP(in_channels=mlp_hidden_channels*len(edge_length_list), hidden_channels=mlp_hidden_channels,
                                                    num_layers=mlp_num_layers, out_channels=mlp_hidden_channels))
            
        self.final_mlp_1 = MLP(in_channels=mlp_hidden_channels*len(edge_lengths), hidden_channels=mlp_hidden_channels,
                               num_layers=mlp_num_layers, out_channels=mlp_hidden_channels)
        self.final_mlp_2 = MLP(in_channels=mlp_hidden_channels, hidden_channels=mlp_hidden_channels,
                               num_layers=mlp_num_layers, out_channels=output_features)
    
    def forward(self, data, plot=False):
    # def forward(self, data):
        # plot=False
        x, batch = data.x, data.batch
        outer_loop_network_outputs = []

        network_idx = 0

        for outer_idx, edge_length_list in enumerate(self.edge_lengths):
            inner_loop_network_outputs = []
            for edge_length in edge_length_list:

                # edge_length_attr_name = f"edge_index_{edge_length}_locations"
                edge_index_locations = abs(data.edge_index[0,:] - data.edge_index[1,:])==edge_length
                edge_index_locations += (data.edge_index[0,:] - data.edge_index[1,:]) == 0
                # edge_index_locations = getattr(data, edge_length_attr_name, None)
                edge_index = data.edge_index[:,edge_index_locations]
                if edge_index is None:
                    raise ValueError(f"Edge length of {edge_length} not found in data")
                
                current_network = self.edge_length_networks[network_idx]

                x_cur = current_network.convs[0](x=x, edge_index=edge_index, edge_weight=None)
                batch_cur = batch

                x_pooling = []
                x_pooling.append(x_cur)

                pooling_indices = []
                pooling_edge_indices = []

                # if plot:
                #     if network_idx%5 == 0:
                #         cur_data = Data(x=x_cur, edge_index=edge_index)
                #         plt.figure(figsize=(14,12))
                #         g = torch_geometric.utils.to_networkx(cur_data, to_undirected=True, remove_self_loops=True)
                #         plotting_labels = {node: str(node) for node in g.nodes()}
                #         plotting_pos = nx.circular_layout(g)
                #         nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                #         plt.title(f"Original, edge length: {edge_length}")
                #         plt.show()

                for i, pooler in enumerate(current_network.poolers):
                    if current_network.pool_types[i] == "ASA":
                        x_cur, edge_index, edge_weight, batch_cur, index = pooler(
                            x=x_cur, edge_index=edge_index, batch=batch_cur)
                    else:
                        x_cur, edge_index, edge_weight, batch_cur, index, score = pooler(
                            x=x_cur, edge_index=edge_index, batch=batch_cur)
                    
                    pooling_indices.append(index)
                    pooling_edge_indices.append(edge_index)
                    
                    if plot:
                        if network_idx%5 == 0:
                            plt.figure(figsize=(14,12))
                            plot_data = Data(x=x_cur, edge_index=edge_index)
                            g = torch_geometric.utils.to_networkx(plot_data, to_undirected=True, remove_self_loops=True)
                            plotting_labels = {node: plotting_labels[index[node].item()] for node in g.nodes()}
                            plotting_pos = {node: plotting_pos[index[node].item()] for node in g.nodes()}
                            nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                            plt.title(f"Pooled level {i}, edge length: {edge_length}")
                            plt.show()
                    
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
            
            x_cur = torch.cat(inner_loop_network_outputs, dim=1)
            outer_loop_network_outputs.append(self.edge_lengths_final_mlps[outer_idx](x_cur))

        x = torch.cat(outer_loop_network_outputs, dim=1)
        x = self.final_mlp_1(x)
        batch = data.batch
        for cluster in self.cluster_list:
            cur_cluster = cluster.repeat(len(batch)//len(cluster)) + (torch.arange(len(batch))//len(cluster))*(torch.max(cluster)+1)
            x, batch = max_pool_x(cur_cluster, x, batch)
        x = global_mean_pool(x, batch)
        x = self.final_mlp_2(x)
        return x
    

class GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling(Module):
    def __init__(self, num_node_features, mp_hidden_channels, mlp_hidden_channels, num_layers, output_features, pool_ratios, pool_types, 
                 mlp_num_layers, edge_lengths, jk:str, mp_layer=GraphSAGE, num_nodes:int=50, final_dropout=0.0):
        super(GraphSAGESeparateEdgeLengthsWithMLPsNoFinalPooling, self).__init__()
        self.edge_lengths = edge_lengths
        self.edge_length_networks = ModuleList()
        self.edge_lengths_final_mlps = ModuleList()
        self.num_nodes = num_nodes

        for edge_length_list in edge_lengths:
            for _ in edge_length_list:
                self.edge_length_networks.append(CustomUNet(num_node_features=num_node_features, 
                                                            mp_hidden_channels=mp_hidden_channels,
                                                            mlp_hidden_channels=mlp_hidden_channels,
                                                            num_layers=num_layers,
                                                            output_features=mlp_hidden_channels,
                                                            pool_ratios=pool_ratios,
                                                            pool_types=pool_types,
                                                            mlp_num_layers=mlp_num_layers,
                                                            jk=jk,
                                                            mp_layer=mp_layer))
            self.edge_lengths_final_mlps.append(MLP(in_channels=mlp_hidden_channels*len(edge_length_list), hidden_channels=mlp_hidden_channels,
                                                    num_layers=mlp_num_layers, out_channels=mlp_hidden_channels))
            
        self.final_mlp_1 = MLP(in_channels=mlp_hidden_channels*len(edge_lengths), hidden_channels=mlp_hidden_channels,
                               num_layers=mlp_num_layers, out_channels=mlp_hidden_channels, dropout=final_dropout)
        self.final_mlp_2 = MLP(in_channels=mlp_hidden_channels*self.num_nodes, hidden_channels=mlp_hidden_channels,
                               num_layers=mlp_num_layers, out_channels=output_features, dropout=final_dropout)
    
    def forward(self, data, plot=False):
    # def forward(self, data):
        # plot=False
        x, batch = data.x, data.batch
        outer_loop_network_outputs = []

        network_idx = 0

        for outer_idx, edge_length_list in enumerate(self.edge_lengths):
            inner_loop_network_outputs = []
            for edge_length in edge_length_list:

                # edge_length_attr_name = f"edge_index_{edge_length}_locations"
                edge_index_locations = abs(data.edge_index[0,:] - data.edge_index[1,:])==edge_length
                edge_index_locations += (data.edge_index[0,:] - data.edge_index[1,:]) == 0
                # edge_index_locations = getattr(data, edge_length_attr_name, None)
                edge_index = data.edge_index[:,edge_index_locations]
                if edge_index is None:
                    raise ValueError(f"Edge length of {edge_length} not found in data")
                
                current_network = self.edge_length_networks[network_idx]

                x_cur = current_network.convs[0](x=x, edge_index=edge_index, edge_weight=None)
                batch_cur = batch

                x_pooling = []
                x_pooling.append(x_cur)

                pooling_indices = []
                pooling_edge_indices = []

                # if plot:
                #     if network_idx%5 == 0:
                #         cur_data = Data(x=x_cur, edge_index=edge_index)
                #         plt.figure(figsize=(14,12))
                #         g = torch_geometric.utils.to_networkx(cur_data, to_undirected=True, remove_self_loops=True)
                #         plotting_labels = {node: str(node) for node in g.nodes()}
                #         plotting_pos = nx.circular_layout(g)
                #         nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                #         plt.title(f"Original, edge length: {edge_length}")
                #         plt.show()

                for i, pooler in enumerate(current_network.poolers):
                    if current_network.pool_types[i] == "ASA":
                        x_cur, edge_index, edge_weight, batch_cur, index = pooler(
                            x=x_cur, edge_index=edge_index, batch=batch_cur)
                    else:
                        x_cur, edge_index, edge_weight, batch_cur, index, score = pooler(
                            x=x_cur, edge_index=edge_index, batch=batch_cur)
                    
                    pooling_indices.append(index)
                    pooling_edge_indices.append(edge_index)
                    
                    if plot:
                        if network_idx%5 == 0:
                            plt.figure(figsize=(14,12))
                            plot_data = Data(x=x_cur, edge_index=edge_index)
                            g = torch_geometric.utils.to_networkx(plot_data, to_undirected=True, remove_self_loops=True)
                            plotting_labels = {node: plotting_labels[index[node].item()] for node in g.nodes()}
                            plotting_pos = {node: plotting_pos[index[node].item()] for node in g.nodes()}
                            nx.draw(g, labels=plotting_labels, pos=plotting_pos)
                            plt.title(f"Pooled level {i}, edge length: {edge_length}")
                            plt.show()
                    
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
            
            x_cur = torch.cat(inner_loop_network_outputs, dim=1)
            outer_loop_network_outputs.append(self.edge_lengths_final_mlps[outer_idx](x_cur))

        x = torch.cat(outer_loop_network_outputs, dim=1)
        x = self.final_mlp_1(x)
        x = x.reshape(torch.max(data.batch)+1, -1)
        x = self.final_mlp_2(x)
        return x