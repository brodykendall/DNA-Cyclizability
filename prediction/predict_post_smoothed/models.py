from imports import *
import matplotlib.pyplot as plt
import networkx as nx
    
class FeedForward(Module):
    def __init__(self, nodes_per_layer, num_layers, output_features, dropout, input_features):
        super(FeedForward, self).__init__()

        self.lins = nn.Sequential()
        for i in range(num_layers):
            if i == 0:
                in_features = input_features
                # in_features = 3
                # in_features = 5
            else:
                in_features = nodes_per_layer
            if i == num_layers-1:
                out_features = output_features
            else:
                out_features = nodes_per_layer
            self.lins.append(nn.Linear(in_features, out_features))
            self.lins.append(Dropout(dropout))
            if i < num_layers-1:
                self.lins.append(nn.ReLU())         
    
    def forward(self, x):
        x = self.lins(x)
        return x