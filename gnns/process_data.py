from helper_functions import *

def get_data_from_csv(filepath):
    dat = pd.read_csv(filepath,delimiter = ",")
    data_list = []
    ys = find_c0new(dat)
    for i, sequence_nt in enumerate(dat["Sequence"]):
        x = dnaOneHot(sequence_nt)
        edge_index = get_forward_and_backward_edges(len(sequence_nt))
        y = ys[i,:].reshape(1, -1)
        data_list.append(Data(x=x, y=y, edge_index=edge_index))
    return data_list

