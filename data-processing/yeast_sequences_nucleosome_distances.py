from Bio import SeqIO
import pandas as pd

# Read the yeast genome from the FASTA file
# genome_file = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/sgdSAC2.fasta"

# Downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/
genome_file = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/sacCer3.fa"
with open(genome_file, "r") as f:
    genome_data = f.read().split("\n")

# Parse the nucleosome positions from Supplementary Table 2
nucleosome_file = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/sgdSAC3_ncp.txt"
nucleosome_positions = pd.read_csv(nucleosome_file, delim_whitespace=True, header=None, names=['chrID', 'position', 'ncp_score', 'ncp_ratio'])

sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Iterate over yeast genome
half_window_size=200
subsequence_length = 50
for chromosome in sequences:
    chromosome_df_list = []
    # if chromosome == "chrV":
    #     continue
    print(chromosome)
    full_chr_sequence = str(sequences[chromosome].seq)
    for nuc_pos in nucleosome_positions.loc[nucleosome_positions['chrID']==chromosome].position:
        center_starting_point = nuc_pos - 25
        for j in range(center_starting_point-half_window_size, center_starting_point+half_window_size):
            subsequence = full_chr_sequence[j:j + subsequence_length]
            distance = j-center_starting_point
            chromosome_df_list.append({"sequence": subsequence, "distance_to_nucleosome": distance, "chrID": chromosome, 
                                       "nucleosome_center": nuc_pos, "position": j})
            
    # for j in range(len(full_chr_sequence) - window_size + 1):
    #     subsequence = full_chr_sequence[j:j + window_size]
    #     distance = closest_nucleosome(chromosome, j)
    #     chromosome_df_list.append({"sequence": subsequence, "distance_to_nucleosome": distance, "chrID": chromosome, 'position': j})

    chromosome_df = pd.DataFrame(chromosome_df_list)

    chromosome_df.to_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/yeast_{chromosome}_nucleosome_alignment_window_size{half_window_size}.csv", index=False)