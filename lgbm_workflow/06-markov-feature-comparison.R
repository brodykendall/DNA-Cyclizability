library(tidyverse)
library(markovchain)
library(modelr)

#Sequence functions
source("scripts/functions/sequence-functions.R")

#Markov feature functions
source("scripts/functions/markov-functions.R")

#Set seed for reproducibility
set.seed(50)

kfold_obj <- readRDS("data/Created/train-10-fold.rds")

likelihood_ratio_outcomes <- list()

for(i in 1:3){
  
  #Calculate likelihood ratio using different quantile cutoffs for
  #low and high cyclizability scores
  kfold_obj_loop <- kfold_obj %>%
    mutate(cutoffs_90_10 = map(train, ~ quantile(.x$C0, c(0.1, 0.9))),
           cutoffs_80_20 = map(train, ~ quantile(.x$C0, c(0.2, 0.8))),
           cutoffs_66_33 = map(train, ~ quantile(.x$C0, c(0.33, 0.67))),
           cutoffs_50_50 = map(train, ~ quantile(.x$C0, c(0.5, 0.5))),
           eff_mat0 = map2(train, cutoffs_90_10,
                           ~ position_specific_score_matrices(.x[.x$C0 > .y[2],], order = i)),
           ineff_mat0 = map2(train, cutoffs_90_10,
                             ~ position_specific_score_matrices(.x[.x$C0 < .y[1],], order = i)),
           eff_mat1 = map2(train, cutoffs_80_20,
                           ~ position_specific_score_matrices(.x[.x$C0 > .y[2],], order = i)),
           ineff_mat1 = map2(train, cutoffs_80_20,
                             ~ position_specific_score_matrices(.x[.x$C0 < .y[1],], order = i)),
           eff_mat2 = map2(train, cutoffs_66_33,
                           ~ position_specific_score_matrices(.x[.x$C0 > .y[2],], order = i)),
           ineff_mat2 = map2(train, cutoffs_66_33,
                             ~ position_specific_score_matrices(.x[.x$C0 < .y[1],], order = i)),
           eff_mat3 = map2(train, cutoffs_50_50,
                           ~ position_specific_score_matrices(.x[.x$C0 > .y[2],], order = i)),
           ineff_mat3 = map2(train, cutoffs_50_50,
                             ~ position_specific_score_matrices(.x[.x$C0 < .y[1],], order = i)),
           ratio_mat0 = map2(eff_mat0, ineff_mat0, ~ map2(.x, .y, `-`)),
           ratio_mat1 = map2(eff_mat1, ineff_mat1, ~ map2(.x, .y, `-`)),
           ratio_mat2 = map2(eff_mat2, ineff_mat2, ~ map2(.x, .y, `-`)),
           ratio_mat3 = map2(eff_mat3, ineff_mat3, ~ map2(.x, .y, `-`)),
           train_ratio_score_0 = map2(train, ratio_mat0, pssm_sum_score),
           train_ratio_score_1 = map2(train, ratio_mat1, pssm_sum_score),
           train_ratio_score_2 = map2(train, ratio_mat2, pssm_sum_score),
           train_ratio_score_3 = map2(train, ratio_mat3, pssm_sum_score),
           test_ratio_score_0 = map2(test, ratio_mat0, pssm_sum_score),
           test_ratio_score_1 = map2(test, ratio_mat1, pssm_sum_score),
           test_ratio_score_2 = map2(test, ratio_mat2, pssm_sum_score),
           test_ratio_score_3 = map2(test, ratio_mat3, pssm_sum_score),
           train = pmap(list(train, train_ratio_score_0, train_ratio_score_1, train_ratio_score_2, train_ratio_score_3),
                        ~ bind_cols(..1, data.frame(ratio_score_0 = ..2,
                                                    ratio_score_1 = ..3,
                                                    ratio_score_2 = ..4,
                                                    ratio_score_3 = ..5))),
           test = pmap(list(test, test_ratio_score_0, test_ratio_score_1, test_ratio_score_2, test_ratio_score_3),
                       ~ bind_cols(..1, data.frame(ratio_score_0 = ..2,
                                                   ratio_score_1 = ..3,
                                                   ratio_score_2 = ..4,
                                                   ratio_score_3 = ..5))))
  
  #Calculate mean Pearson correlation across folds
  ratio_cutoff_comparison <- kfold_obj_loop %>%
    mutate(markov_ratio_pearson = map(test,
                                      ~ data.frame(r0_pc = cor(.x$C0, .x$ratio_score_0, method = "pearson"),
                                                   r1_pc = cor(.x$C0, .x$ratio_score_1, method = "pearson"),
                                                   r2_pc = cor(.x$C0, .x$ratio_score_2, method = "pearson"),
                                                   r3_pc = cor(.x$C0, .x$ratio_score_3, method = "pearson")))) %>%
    pull(markov_ratio_pearson) %>%
    bind_rows() %>%
    summarize(avg0 = mean(r0_pc),
              avg1 = mean(r1_pc),
              avg2 = mean(r2_pc),
              avg3 = mean(r3_pc),
              sd0 = sd(r0_pc),
              sd1 = sd(r1_pc),
              sd2 = sd(r2_pc),
              sd3 = sd(r3_pc))
  
  #Gather results
  ratio_cutoff_comparison <- bind_cols(cutoff = c("q90/10", "q80/20", "q67/33", "q50/50"),
                                       gather(ratio_cutoff_comparison[,1:4], value = "avg_pearson"),
                                       gather(ratio_cutoff_comparison[,5:8], value = "sd")) %>%
    select(-`key...2`, -`key...4`) %>%
    mutate(order = i)
  
  likelihood_ratio_outcomes[[i]] <- ratio_cutoff_comparison
  
}

likelihood_ratio_outcomes_combined<- likelihood_ratio_outcomes %>%
  bind_rows()

head(likelihood_ratio_outcomes_combined)

likelihood_ratio_outcomes_combined_display <- likelihood_ratio_outcomes_combined %>%
  select(order, cutoff, avg_pearson, sd) %>%
  arrange(order, cutoff) %>%
  mutate(avg_pearson = round(avg_pearson, 4),
         sd = round(sd, 4))

write_csv(likelihood_ratio_outcomes_combined, "results/ratio-comparisons.csv")
write_csv(likelihood_ratio_outcomes_combined_display, "results/ratio-comparisons-display.csv")

#Add optimal ratio to kfold obj
kfold_obj <- kfold_obj %>%
  mutate(cutoffs_80_20 = map(train, ~ quantile(.x$C0, c(0.2, 0.8))),
         eff_mat2 = map2(train, cutoffs_80_20,
                         ~ position_specific_score_matrices(.x[.x$C0 > .y[2],], order = 2)),
         ineff_mat2 = map2(train, cutoffs_80_20,
                           ~ position_specific_score_matrices(.x[.x$C0 < .y[1],], order = 2)),
         ratio_mat2 = map2(eff_mat2, ineff_mat2, ~ map2(.x, .y, `-`)),
         train_ratio_score_2 = map2(train, ratio_mat2, pssm_sum_score),
         test_ratio_score_2 = map2(test, ratio_mat2, pssm_sum_score),
         train = pmap(list(train, train_ratio_score_2),
                      ~ bind_cols(..1, data.frame(ratio_score_2nd = ..2))),
         test = pmap(list(test, test_ratio_score_2),
                     ~ bind_cols(..1, data.frame(ratio_score_2nd = ..2))))

#Save kfold object with ratios
saveRDS(kfold_obj, "data/Created/train-10-fold-ratio.rds")
