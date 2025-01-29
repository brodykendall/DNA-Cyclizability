from imports import *
from torch.utils.data import Dataset
from helper_functions import *

class DNASeqDataset(Dataset):
    def __init__(self, src_file, edge_lengths, kernel_length):
        unprocessed_data = pd.read_csv(src_file,delimiter = ",")
        base_x = []
        for sequence_nt in unprocessed_data["Sequence"]:
            base_x.append(dnaOneHot(sequence_nt).reshape(1,50,4))
        base_x = torch.cat(base_x)
        base_x = base_x.reshape((base_x.shape[0],50,4,1))
        # base_x = np.transpose(base_x, (0, 3, 1, 2))
        base_x = base_x.permute(0, 3, 1, 2)
        self.x = []
        for edge_length_list in edge_lengths:
            for edge_length in edge_length_list:
                x_cur = [torch.zeros_like(base_x[:,:,range(kernel_length)])]
                for i in range(edge_length):
                    x_cur.append(base_x[:,:,range(i,50,edge_length)])
                    x_cur.append(torch.zeros_like(base_x[:,:,range(kernel_length)]))
                x_cur = torch.cat(x_cur, dim=2)
                self.x.append(x_cur)

        self.y = find_c0new(unprocessed_data)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return [edge_length_tensor[idx] for edge_length_tensor in self.x], self.y[idx]