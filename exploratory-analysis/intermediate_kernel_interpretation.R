
examine_3x4kernels = function(row) {
  ret = c()
  for(i in 1:4){
    for(j in 1:4){
      for(k in 1:4){
        val = row[i] + row[4+j] + row[8+k] + row[13]
        if(val > 0) {
          ret = append(ret, paste0("nonzero: ", nucleotides[i], ", ", nucleotides[j], ", ", nucleotides[k], ": ", val))
        }
      }
    }
  }
  return(ret)
}

# 16 3x4 kernels:
tiling_first_conv_kernels.conv_only_2 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_2.rds")
tiling_first_conv_biases.conv_only_2 = readRDS("data/Created/tiling_first_conv_biases_conv_only_2.rds")
tiling_first_conv_kernels_full.conv_only_2 = cbind(tiling_first_conv_kernels.conv_only_2, tiling_first_conv_biases.conv_only_2)
apply(tiling_first_conv_kernels_full.conv_only_2, 1, examine_3x4kernels)

# 4 3x4 kernels:
tiling_first_conv_kernels.conv_only_4 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_4.rds")
tiling_first_conv_biases.conv_only_4 = readRDS("data/Created/tiling_first_conv_biases_conv_only_4.rds")
tiling_first_conv_kernels_full.conv_only_4 = cbind(tiling_first_conv_kernels.conv_only_4, tiling_first_conv_biases.conv_only_4)
apply(tiling_first_conv_kernels_full.conv_only_4, 1, examine_3x4kernels)

# 8 3x4 kernels:
tiling_first_conv_kernels.conv_only_5 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_5.rds")
tiling_first_conv_biases.conv_only_5 = readRDS("data/Created/tiling_first_conv_biases_conv_only_5.rds")
tiling_first_conv_kernels_full.conv_only_5 = cbind(tiling_first_conv_kernels.conv_only_5, tiling_first_conv_biases.conv_only_5)
apply(tiling_first_conv_kernels_full.conv_only_5, 1, examine_3x4kernels)


examine_5x4kernels = function(row) {
  ret = c()
  for(i in 1:4){
    for(j in 1:4){
      for(k in 1:4){
        for(l in 1:4){
          for(m in 1:4){
            val = row[i] + row[4+j] + row[8+k] + row[12+l] + row[16+m] + row[17]
            if(val > 0) {
              ret = append(ret, paste0("nonzero: ", nucleotides[i], ", ", nucleotides[j], ", ", 
                                       nucleotides[k], ", ", nucleotides[l], ", ", nucleotides[m], ": ", val))
            }
          }
        }
      }
    }
  }
  return(ret)
}


# 48 5x4 kernels:
tiling_first_conv_kernels.conv_only_3 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_3.rds")
tiling_first_conv_biases.conv_only_3 = readRDS("data/Created/tiling_first_conv_biases_conv_only_3.rds")
tiling_first_conv_kernels_full.conv_only_3 = cbind(tiling_first_conv_kernels.conv_only_3, tiling_first_conv_biases.conv_only_3)
apply(tiling_first_conv_kernels_full.conv_only_3, 1, examine_5x4kernels)$feature19









#### Dinucleotide input:

dinucleotide_ordered = c("AA", "TT", "AT", "TA",
                         "AC", "GT", "TC", "GA",
                         "AG", "CT", "TG", "CA",
                         "CG", "GC", "CC", "GG")

# 8 21x16 kernels:
tiling_first_conv_kernels.conv_only_di_5 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_di_5.rds")
tiling_first_conv_biases.conv_only_di_5 = readRDS("data/Created/tiling_first_conv_biases_conv_only_di_5.rds")
tiling_first_conv_kernels_full.conv_only_di_5 = cbind(tiling_first_conv_kernels.conv_only_di_5, tiling_first_conv_biases.conv_only_di_5)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_di_5.formatted_features_list = 1:8 %>%
  map(~expand.grid(dinucleotide=dinucleotides, position=1:21)) %>%
  map2(1:8, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_di_5[.y,])))) %>%
  map(~.x %>%
        mutate(across(dinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = dinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_di_5.heatmaps = tiling_first_conv_kernels.conv_only_di_5.formatted_features_list %>%
  map2(1:8, ~ggplot(.x, aes(position, dinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_di_5.heatmaps[[4]]

tiling_first_conv_kernels.conv_only_di_5.heatmaps[[8]]

tiling_first_conv_kernels.conv_only_di_5.heatmaps[[5]]

tiling_first_conv_kernels.conv_only_di_5.heatmaps[[6]]






# 48 21x16 kernels:
tiling_first_conv_kernels.conv_only_di = readRDS("data/Created/tiling_first_conv_kernels_conv_only_di.rds")
tiling_first_conv_biases.conv_only_di = readRDS("data/Created/tiling_first_conv_biases_conv_only_di.rds")
tiling_first_conv_kernels_full.conv_only_di = cbind(tiling_first_conv_kernels.conv_only_di, tiling_first_conv_biases.conv_only_di)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_di.formatted_features_list = 1:48 %>%
  map(~expand.grid(dinucleotide=dinucleotides, position=1:21)) %>%
  map2(1:48, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_di[.y,])))) %>%
  map(~.x %>%
        mutate(across(dinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = dinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_di.heatmaps = tiling_first_conv_kernels.conv_only_di.formatted_features_list %>%
  map2(1:48, ~ggplot(.x, aes(position, dinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_di.heatmaps[[1]]

tiling_first_conv_kernels.conv_only_di.heatmaps[[48]]









#### Trinucleotide input:




# 4 21x64 kernels:
tiling_first_conv_kernels.conv_only_tri_4 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_tri_4.rds")
tiling_first_conv_biases.conv_only_tri_4 = readRDS("data/Created/tiling_first_conv_biases_conv_only_tri_4.rds")
tiling_first_conv_kernels_full.conv_only_tri_4 = cbind(tiling_first_conv_kernels.conv_only_tri_4, tiling_first_conv_biases.conv_only_tri_4)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_tri_4.formatted_features_list = 1:4 %>%
  map(~expand.grid(trinucleotide=trinucleotides, position=1:21)) %>%
  map2(1:4, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_tri_4[.y,])))) %>%
  map(~.x %>%
        mutate(across(trinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = trinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_tri_4.heatmaps = tiling_first_conv_kernels.conv_only_tri_4.formatted_features_list %>%
  map2(1:4, ~ggplot(.x, aes(position, trinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_tri_4.heatmaps[[3]]

tiling_first_conv_kernels.conv_only_tri_4.heatmaps[[1]]

tiling_first_conv_kernels.conv_only_tri_4.heatmaps[[2]]

tiling_first_conv_kernels.conv_only_tri_4.heatmaps[[4]]










# 8 21x64 kernels:
tiling_first_conv_kernels.conv_only_tri_5 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_tri_5.rds")
tiling_first_conv_biases.conv_only_tri_5 = readRDS("data/Created/tiling_first_conv_biases_conv_only_tri_5.rds")
tiling_first_conv_kernels_full.conv_only_tri_5 = cbind(tiling_first_conv_kernels.conv_only_tri_5, tiling_first_conv_biases.conv_only_tri_5)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_tri_5.formatted_features_list = 1:8 %>%
  map(~expand.grid(trinucleotide=trinucleotides, position=1:21)) %>%
  map2(1:8, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_tri_5[.y,])))) %>%
  map(~.x %>%
        mutate(across(trinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = trinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_tri_5.heatmaps = tiling_first_conv_kernels.conv_only_tri_5.formatted_features_list %>%
  map2(1:8, ~ggplot(.x, aes(position, trinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[5]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[6]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[1]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[2]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[3]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[4]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[7]]

tiling_first_conv_kernels.conv_only_tri_5.heatmaps[[8]]









# 48 21x64 kernels:
tiling_first_conv_kernels.conv_only_tri = readRDS("data/Created/tiling_first_conv_kernels_conv_only_tri.rds")
tiling_first_conv_biases.conv_only_tri = readRDS("data/Created/tiling_first_conv_biases_conv_only_tri.rds")
tiling_first_conv_kernels_full.conv_only_tri = cbind(tiling_first_conv_kernels.conv_only_tri, tiling_first_conv_biases.conv_only_tri)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_tri.formatted_features_list = 1:48 %>%
  map(~expand.grid(trinucleotide=trinucleotides, position=1:21)) %>%
  map2(1:48, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_tri[.y,])))) %>%
  map(~.x %>%
        mutate(across(trinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = trinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_tri.heatmaps = tiling_first_conv_kernels.conv_only_tri.formatted_features_list %>%
  map2(1:48, ~ggplot(.x, aes(position, trinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_tri.heatmaps[[19]]

tiling_first_conv_kernels.conv_only_tri.heatmaps[[31]]

tiling_first_conv_kernels.conv_only_tri.heatmaps[[25]]







# 48 31x64 kernels:
tiling_first_conv_kernels.conv_only_tri_2 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_tri_2.rds")
tiling_first_conv_biases.conv_only_tri_2 = readRDS("data/Created/tiling_first_conv_biases_conv_only_tri_2.rds")
tiling_first_conv_kernels_full.conv_only_tri_2 = cbind(tiling_first_conv_kernels.conv_only_tri_2, tiling_first_conv_biases.conv_only_tri_2)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_tri_2.formatted_features_list = 1:48 %>%
  map(~expand.grid(trinucleotide=trinucleotides, position=1:31)) %>%
  map2(1:48, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_tri_2[.y,])))) %>%
  map(~.x %>%
        mutate(across(trinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = trinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_tri_2.heatmaps = tiling_first_conv_kernels.conv_only_tri_2.formatted_features_list %>%
  map2(1:48, ~ggplot(.x, aes(position, trinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_tri_2.heatmaps[[8]]

tiling_first_conv_kernels.conv_only_tri_2.heatmaps[[15]]

tiling_first_conv_kernels.conv_only_tri_2.heatmaps[[2]]

tiling_first_conv_kernels.conv_only_tri_2.heatmaps[[9]]





trinucleotide_ordered = c("AAA", "TTT", "AAT", "ATT",
                          "ATA", "TAT", "TAA", "TTA",
                          "TAC", "GTA",
                          "GAA", "TTC", "CAG", "CTG",
                          "AAC", "GTT", "GGA", "TCC",
                          "AAG", "CTT", "TAG", "CTA",
                          "TCT", "AGA", 
                          "CAA", "TTG", "CAT", "ATG",
                          
                          "GAT", "ATC",
                          "CAC", "GTG", "CTC", "GAG",
                          "GAC", "GTC",
                          "TCA", "TGA", "AGT", "ACT",
                          "TGT", "ACA",
                          "CGA", "TCG", 
                          "CGT", "ACG", "GGT", "ACC",
                          "AGC", "GCT", "AGG", "CCT",
                          "TGC", "GCA", "TGG", "CCA",
                          "GGC", "GCC", "GCG", "CGC",
                          "CGG", "CCG", "CCC", "GGG")

# 48 41x64 kernels:
tiling_first_conv_kernels.conv_only_tri_3 = readRDS("data/Created/tiling_first_conv_kernels_conv_only_tri_3.rds")
tiling_first_conv_biases.conv_only_tri_3 = readRDS("data/Created/tiling_first_conv_biases_conv_only_tri_3.rds")
tiling_first_conv_kernels_full.conv_only_tri_3 = cbind(tiling_first_conv_kernels.conv_only_tri_3, tiling_first_conv_biases.conv_only_tri_3)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_tri_3.formatted_features_list = 1:48 %>%
  map(~expand.grid(trinucleotide=trinucleotides, position=1:41)) %>%
  map2(1:48, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_tri_3[.y,])))) %>%
  map(~.x %>%
        mutate(across(trinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = trinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_tri_3.heatmaps = tiling_first_conv_kernels.conv_only_tri_3.formatted_features_list %>%
  map2(1:48, ~ggplot(.x, aes(position, trinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_tri_3.heatmaps[[48]]

tiling_first_conv_kernels.conv_only_tri_3.heatmaps[[6]]

tiling_first_conv_kernels.conv_only_tri_3.heatmaps[[45]]








#### Prediction on Residuals of model trained on Fourier features
# 8 21x64 kernels:
tiling_first_conv_kernels.conv_only_tri_5_fourier_resid = readRDS("data/Created/tiling_first_conv_kernels_conv_only_tri_5_fourier_resid.rds")
tiling_first_conv_biases.conv_only_tri_5_fourier_resid = readRDS("data/Created/tiling_first_conv_biases_conv_only_tri_5_fourier_resid.rds")
tiling_first_conv_kernels_full.conv_only_tri_5_fourier_resid = cbind(tiling_first_conv_kernels.conv_only_tri_5_fourier_resid, tiling_first_conv_biases.conv_only_tri_5_fourier_resid)

# Feature heatmaps:
tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.formatted_features_list = 1:8 %>%
  map(~expand.grid(trinucleotide=trinucleotides, position=1:21)) %>%
  map2(1:8, ~.x %>%
         mutate(weight = as.vector(t(tiling_first_conv_kernels.conv_only_tri_5_fourier_resid[.y,])))) %>%
  map(~.x %>%
        mutate(across(trinucleotide,
                      ~factor(.x, ordered = TRUE,
                              levels = trinucleotide_ordered))))

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps = tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.formatted_features_list %>%
  map2(1:8, ~ggplot(.x, aes(position, trinucleotide)) + 
         geom_tile(aes(fill=weight)) +
         scale_fill_gradient2("Positive Weight", low="red", mid="white", high="green", midpoint=0) +
         ggtitle(paste0("Kernel ", .y)))

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[5]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[6]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[1]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[2]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[3]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[4]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[7]]

tiling_first_conv_kernels.conv_only_tri_5_fourier_resid.heatmaps[[8]]
