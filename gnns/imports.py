import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch_geometric
import torch_geometric.nn as tg_nn

from torch.optim.lr_scheduler import ReduceLROnPlateau

from torch.nn import Module
from torch.nn import ModuleList

from torch_geometric.nn import MLP
from torch_geometric.nn import GraphSAGE
from torch_geometric.nn import GCN
from torch_geometric.nn import GAT
from torch_geometric.nn import SAGPooling
from torch_geometric.nn import TopKPooling
from torch_geometric.nn import ASAPooling
from torch_geometric.nn import EdgePooling
from torch_geometric.nn import global_mean_pool
from torch_geometric.nn import global_add_pool
from torch_geometric.nn import avg_pool
from torch_geometric.nn import max_pool
from torch_geometric.nn import avg_pool_x
from torch_geometric.nn import max_pool_x

from torch_geometric.data import Data
from torch_geometric.data import InMemoryDataset
from torch_geometric.loader import DataLoader