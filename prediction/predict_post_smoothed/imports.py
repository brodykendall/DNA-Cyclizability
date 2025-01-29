import numpy as np
import pandas as pd

import torch
import torch.nn as nn

from torch.optim.lr_scheduler import ReduceLROnPlateau

from torch.nn import Module
from torch.nn import ModuleList

from torch.nn import Conv2d
from torch.nn import Linear
from torch.nn import ReLU
# from torch.nn import Dropout1d
# from torch.nn import Dropout2d
from torch.nn import Dropout

from torch.utils.data import DataLoader