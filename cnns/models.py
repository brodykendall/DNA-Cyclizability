from imports import *
import matplotlib.pyplot as plt
import networkx as nx
    
class ConvSeparateEdgeLengthsWithMLPs(Module):
    def __init__(self, num_node_features, conv_hidden_channels, mlp_hidden_channels, output_features, 
                 edge_lengths, kernel_length, conv_dropout, mlp_dropout):
        super(ConvSeparateEdgeLengthsWithMLPs, self).__init__()
        self.edge_lengths = edge_lengths
        self.edge_length_networks = ModuleList()
        self.edge_lengths_final_mlps = ModuleList()

        for edge_length_list in edge_lengths:
            cur_outputs = 0
            for edge_length in edge_length_list:
                self.edge_length_networks.append(nn.Sequential(
                    Conv2d(in_channels=1, out_channels=conv_hidden_channels, kernel_size=(kernel_length,num_node_features)),
                    Dropout(conv_dropout)))
                cur_outputs += (51-kernel_length)+(kernel_length)*(edge_length+1)
            self.edge_lengths_final_mlps.append(nn.Sequential(
                Linear(in_features=conv_hidden_channels*cur_outputs, out_features=mlp_hidden_channels*len(edge_length_list)),
                ReLU(),
                Linear(in_features=mlp_hidden_channels*len(edge_length_list), out_features=mlp_hidden_channels*len(edge_length_list)),
                ReLU(),
                Linear(in_features=mlp_hidden_channels*len(edge_length_list), out_features=mlp_hidden_channels),
                ReLU(),
                Dropout(mlp_dropout)))
            
        self.final_mlp = nn.Sequential(
            Linear(in_features=mlp_hidden_channels*len(edge_lengths), out_features=output_features))
        # self.final_mlp_2 = nn.Sequential(
        #     Linear(in_features=mlp_hidden_channels, out_features=mlp_hidden_channels),
        #     ReLU(),
        #     Dropout(mlp_dropout),
        #     Linear(in_features=mlp_hidden_channels, out_features=output_features))            
    
    def forward(self, x, plot=False):
        outer_loop_network_outputs = []

        network_idx = 0

        for outer_idx, edge_length_list in enumerate(self.edge_lengths):
            inner_loop_network_outputs = []
            for inner_idx, edge_length in enumerate(edge_length_list):
                
                current_network = self.edge_length_networks[network_idx]

                inner_loop_network_outputs.append(current_network(x[network_idx]))

                # x_cur = []
                # for i in range(edge_length):
                #     x_cur.append(current_network(x[network_idx]))
                
                # x_cur = torch.cat(x_cur, dim=2)

                # inner_loop_network_outputs.append(x_cur)

                network_idx += 1
            
            x_cur = torch.cat(inner_loop_network_outputs, dim=2)
            x_cur = x_cur.reshape(x_cur.shape[0], -1)
            outer_loop_network_outputs.append(self.edge_lengths_final_mlps[outer_idx](x_cur))

        x = torch.cat(outer_loop_network_outputs, dim=1)
        x = self.final_mlp(x)
        # x = self.final_mlp_2(x)
        return x