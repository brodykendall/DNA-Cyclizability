from Bio import SeqIO
import pandas as pd

# Read the yeast genome from the FASTA file
# genome_file = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/sgdSAC2.fasta"

# Downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/
genome_file = "/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/sacCer3.fa"
with open(genome_file, "r") as f:
    genome_data = f.read().split("\n")

sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Iterate over yeast genome
subsequence_length = 50
for chromosome in sequences:
    chromosome_df_list = []
    # if chromosome == "chrV":
    #     continue
    print(chromosome)
    full_chr_sequence = str(sequences[chromosome].seq)
    for j in range(len(full_chr_sequence)-subsequence_length+1):
        subsequence = full_chr_sequence[j:j + subsequence_length]
        chromosome_df_list.append({"sequence": subsequence, "chrID": chromosome, "position": j})

    chromosome_df = pd.DataFrame(chromosome_df_list)

    chromosome_df.to_csv(f"/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/data/Created/yeast_{chromosome}_1bpresolution_subsequence{subsequence_length}.csv", index=False)