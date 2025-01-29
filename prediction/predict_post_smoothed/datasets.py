from imports import *
from torch.utils.data import Dataset
from helper_functions import *

class PostSmoothedDataset(Dataset):
    def __init__(self, x, y):
        self.X = torch.Tensor(x)
        self.y = torch.Tensor(y)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx,:], self.y[idx]