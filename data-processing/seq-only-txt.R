# Nucleosome:
cycle1 <- read_csv("cycle1.txt")

cycle1 = cycle1 %>% select(Sequence)

write.table(cycle1, "cycle1_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Nucleosome (Rotated):
nuc_rotated <- read_csv("data/Created/nuc_rotated.csv")

nuc_rotated = nuc_rotated %>% select(Sequence)

write.table(nuc_rotated, "data/Created/nuc_rotated_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)



# Random:
cycle3 <- read_csv("cycle3.txt")

cycle3 = cycle3 %>% select(Sequence)

write.table(cycle3, "cycle3_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Random (Rotated):
random_rotated <- read_csv("data/Created/random_rotated.csv")

random_rotated = random_rotated %>% select(Sequence)

write.table(random_rotated, "data/Created/random_rotated_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Tiling:
cycle5 <- read_csv("cycle5.txt")

cycle5 = cycle5 %>% select(Sequence)

write.table(cycle5, "cycle5_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Tiling (Rotated):
tiling_rotated <- read_csv("data/Created/tiling_rotated.csv")

tiling_rotated = tiling_rotated %>% select(Sequence)

write.table(tiling_rotated, "data/Created/tiling_rotated_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# ChrV:
cycle6 <- read_csv("cycle6.txt")

cycle6 = cycle6 %>% select(Sequence)

write.table(cycle6, "cycle6_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# ChrV (Rotated):
chrv_rotated <- read_csv("data/Created/chrv_rotated.csv")

chrv_rotated = chrv_rotated %>% select(Sequence)

write.table(chrv_rotated, "data/Created/chrv_rotated_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)



# ChrV 1bp:
chrv_1bp <- read_csv("data/Created/yeast_chrV_1bpresolution_subsequence50.csv")

chrv_1bp = chrv_1bp %>% select(sequence)

write.table(chrv_1bp, "chrv_1bp_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# ChrV nucleosome aligned:
chrv_aligned <- read_csv("data/Created/yeast_chrV_nucleosome_alignment_window_size200.csv")

chrv_aligned = chrv_aligned %>% select(sequence)

write.table(chrv_aligned, "chrv_aligned_seq_only.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
