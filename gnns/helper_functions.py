from imports import *

def dnaOneHot(sequence, include_indices=False):
    seq_array = np.array(list(sequence))
    code = {"A": [0], "C": [1], "G": [2], "T": [3], "N": [4],
            "a": [0], "c": [1], "g": [2], "t": [3], "n": [4]}
    onehot_encoded_seq = []
    for idx, char in enumerate(seq_array):
        if include_indices:
            onehot_encoded = list(np.zeros(5+len(sequence)))
            onehot_encoded[code[char][0]] = 1
            onehot_encoded[idx+5] = 1
            onehot_encoded.pop(4)
            onehot_encoded_seq.append(onehot_encoded)
        else:
            onehot_encoded = list(np.zeros(5))
            onehot_encoded[code[char][0]] = 1
            onehot_encoded_seq.append(onehot_encoded[0:4])

    if include_indices:
        return torch.Tensor(onehot_encoded_seq).reshape(-1, 4+len(sequence))
    else:
        return torch.Tensor(onehot_encoded_seq).reshape(-1, 4)

def get_forward_and_backward_edges(sequence_length, edge_distances:list[int], self_loops:bool=True):
    ret = torch.Tensor().reshape(2,0)
    for e_d in edge_distances:
        ret = torch.cat((ret, 
                         torch.cat((torch.cat((torch.Tensor(list(range(sequence_length-e_d))).reshape(1,-1), torch.tensor(list(range(e_d,sequence_length))).reshape(1,-1))),
                                    torch.cat((torch.tensor(list(range(e_d,sequence_length))).reshape(1,-1), (torch.Tensor(list(range(sequence_length-e_d))).reshape(1,-1))))), dim=1)), dim=1)
    if self_loops:
        ret = torch.cat((ret, torch.cat((torch.Tensor(list(range(sequence_length))).reshape(1,-1), 
                                         torch.Tensor(list(range(sequence_length))).reshape(1,-1)))), dim=1)
    return ret.type(torch.int64)

def find_c0new(dat):
  mat = np.empty((3,3), float)
  k = 2*np.pi/10.4
  n = np.array([26, 29, 31])
  mat[0:3,0] = 1
  mat[0:3, 1] = np.sin(n*k)
  mat[0:3, 2] = np.cos(n*k)
  inv_mat = np.linalg.inv(mat)
  c0A1A2 = np.array(np.matmul(dat[["n=26", "n=29", "n=31"]], np.transpose(inv_mat)))
  c0Aphi = c0A1A2
  c0Aphi[:,0] = c0A1A2[:,0]
  c0Aphi[:,1] = np.sqrt(c0A1A2[:,1]**2 + c0A1A2[:,2]**2)
  c0Aphi[:,2] <- np.sign(c0A1A2[:,2]) * np.arccos(c0A1A2[:,1]/c0Aphi[:,1])
  ret = torch.Tensor(c0Aphi[:,0]).reshape(-1, 1)
  return ret

def model_improved(cur_loss, best_loss, no_improvement_count, net, model_save_path):
    if cur_loss < best_loss:
        torch.save(net.state_dict(), model_save_path)
        improved = True
        best_loss = cur_loss
        no_improvement_count = 0
    else:
        no_improvement_count += 1
        improved = False
    return improved, best_loss, no_improvement_count