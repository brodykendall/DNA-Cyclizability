from imports import *
from helper_functions import *

class DNASeqGraph(InMemoryDataset):
    def __init__(self, root, raw_graph_name, processed_graph_name, edge_distances:list[int]=[1], 
                 include_indices:bool=False):
        self.raw_graph_name, self.processed_graph_name = raw_graph_name, processed_graph_name
        self.edge_distances = edge_distances
        self.include_indices=include_indices
        super().__init__(root)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        # return ["/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/" + self.raw_graph_name + ".txt"]
        return ["~/Documents/DNA_Cyclizability/" + self.raw_graph_name + ".txt"]

    @property
    def processed_file_names(self):
        return self.processed_graph_name + '.pt'

    def process(self):
        # Read data into huge `Data` list.
        data_list = []

        edge_weight = torch.Tensor()
        for edge_distance in self.edge_distances:
            edge_weight = torch.cat((edge_weight, torch.Tensor([1/(1+edge_distance)]*((50-edge_distance)*2))))
        edge_weight = torch.cat((edge_weight, torch.Tensor([1]*50)))

        for filepath in self.raw_file_names:
            dat = pd.read_csv(filepath,delimiter = ",")
            data_list = []
            ys = find_c0new(dat)
            for i, sequence_nt in enumerate(dat["Sequence"]):
                x = dnaOneHot(sequence_nt, include_indices=self.include_indices)
                edge_index = get_forward_and_backward_edges(len(sequence_nt), self.edge_distances)
                y = ys[i,:].reshape(1, -1)
                data_list.append(Data(x=x, y=y, edge_index=edge_index, edge_weight=edge_weight))

        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])

class DNASeqGraphSeparateEdges(InMemoryDataset):
    def __init__(self, root, raw_graph_dir, raw_graph_name, raw_graph_filetype, variable_of_interest,
                 processed_graph_name, edge_distances:list[int]=[1], include_indices:bool=False):
        self.raw_graph_dir = raw_graph_dir
        self.raw_graph_name, self.processed_graph_name = raw_graph_name, processed_graph_name
        self.raw_graph_file_type = raw_graph_filetype
        self.variable_of_interest = variable_of_interest
        self.edge_distances = edge_distances
        self.include_indices = include_indices
        super().__init__(root)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        # return ["/Users/Brody1/Documents/Northwestern/DNA_Cyclizability/" + self.raw_graph_name + ".txt"]
        # return ["~/Documents/DNA_Cyclizability/" + self.raw_graph_name + ".txt"]
        return [self.raw_graph_dir + self.raw_graph_name + "." + self.raw_graph_file_type]

    @property
    def processed_file_names(self):
        return self.processed_graph_name + '.pt'

    def process(self):
        # Read data into huge `Data` list.
        data_list = []

        for filepath in self.raw_file_names:
            print(f"Processing {filepath}")
            dat = pd.read_csv(filepath,delimiter = ",")
            data_list = []
            # ys = find_c0new(dat)
            ys = torch.Tensor(dat[self.variable_of_interest])
            for i, sequence_nt in enumerate(dat["Sequence"]):
                # Skip first and last sequences:
                if i==0 or i==dat.shape[0]-1:
                    continue
                if i % 1000 == 999:
                    print(f"\t Begin Sequence {i+1}")
                x = dnaOneHot(sequence_nt, include_indices=self.include_indices)
                # y = ys[i,:].reshape(1, -1)
                y = ys[i].reshape(1, -1)
                data_obj = Data(x=x, y=y, edge_index=torch.Tensor())
                for edge_distance in self.edge_distances:
                    edge_index = get_forward_and_backward_edges(len(sequence_nt), [edge_distance], self_loops=False)
                    if data_obj.edge_index.size()[0] != 0:
                        cur_edge_index_locations = torch.arange(data_obj.edge_index.shape[1], 
                                                                data_obj.edge_index.shape[1] + edge_index.shape[1])
                        data_obj.edge_index = torch.cat((data_obj.edge_index, edge_index), dim=1)
                    else:
                        cur_edge_index_locations = torch.arange(edge_index.shape[1])
                        data_obj.edge_index = edge_index
                    setattr(data_obj, f"edge_index_{edge_distance}_locations", cur_edge_index_locations)
                    # setattr(data_obj, f"edge_index_{edge_distance}", edge_index)
                data_list.append(data_obj)

        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])