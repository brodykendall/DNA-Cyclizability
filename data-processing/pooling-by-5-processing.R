tiling_nrow = 82368
set.seed(50)
tiling_train_indices = sample(1:tiling_nrow, tiling_nrow*.9, replace=FALSE)

cycle5 <- read_csv("cycle5.txt")
colnames(cycle5)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")

random_nrow = 12472
set.seed(50)
random_train_indices = sample(1:random_nrow, random_nrow*.9, replace=FALSE)

cycle3 <- read_csv("cycle3.txt")
colnames(cycle3)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

find_pooling_by_5_scores = function(seq, substrs_of_interest=c("A", "T"), offset=0){
  scores = c()
  if (offset != 0) {
    cur_region = substr(seq, 1, offset)
    scores = c(scores, 2*sum(substrs_of_interest %>% 
                               map(~str_count(cur_region, .x)) %>%
                               unlist()) - nchar(cur_region))
  }
  regions = 1:((nchar(seq))%/%5) %>% 
    map(~substr(seq, ((.x-1)*5+offset+1), (.x*5+offset))) %>%
    unlist()
  scores = c(scores, lapply(regions, function(region){
    2*sum(substrs_of_interest %>%
            map(~str_count(region, .x)) %>%
            unlist()) - nchar(region)
  }) %>% unlist())
  return(scores)
}


cycle5_train = cycle5[tiling_train_indices,] %>%
  select(Sequence)
cycle5_test = cycle5[-tiling_train_indices,] %>%
  select(Sequence)

cycle5_pooling_by_5_scores_offset0 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=0)})
cycle5_pooling_by_5_scores_offset1 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=1)})
cycle5_pooling_by_5_scores_offset2 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=2)})
cycle5_pooling_by_5_scores_offset3 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=3)})
cycle5_pooling_by_5_scores_offset4 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=4)})

cycle5_train_pooling_by_5_scores_list = list(cycle5_pooling_by_5_scores_offset0,
                                             cycle5_pooling_by_5_scores_offset1,
                                             cycle5_pooling_by_5_scores_offset2,
                                             cycle5_pooling_by_5_scores_offset3,
                                             cycle5_pooling_by_5_scores_offset4)

cycle5_test_pooling_by_5_scores_offset0 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=0)})
cycle5_test_pooling_by_5_scores_offset1 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=1)})
cycle5_test_pooling_by_5_scores_offset2 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=2)})
cycle5_test_pooling_by_5_scores_offset3 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=3)})
cycle5_test_pooling_by_5_scores_offset4 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=4)})

cycle5_test_pooling_by_5_scores_list = list(cycle5_test_pooling_by_5_scores_offset0,
                                            cycle5_test_pooling_by_5_scores_offset1,
                                            cycle5_test_pooling_by_5_scores_offset2,
                                            cycle5_test_pooling_by_5_scores_offset3,
                                            cycle5_test_pooling_by_5_scores_offset4)

find_alternating_score = function(scores) {
  odd_scores = scores[c(TRUE, FALSE)]
  even_scores = scores[c(FALSE, TRUE)]
  return(abs(sum(odd_scores) - sum(even_scores)))
}

cycle5_train_global_alternating_scores_list = lapply(cycle5_train_pooling_by_5_scores_list,
                                                     function(scores_list) {
                                                       lapply(scores_list, find_alternating_score)
                                                     })
cycle5_train_global_alternating_scores_df = as.data.frame(
  do.call(cbind, cycle5_train_global_alternating_scores_list))

cycle5_test_global_alternating_scores_list = lapply(cycle5_test_pooling_by_5_scores_list,
                                                    function(scores_list) {
                                                      lapply(scores_list, find_alternating_score)
                                                    })
cycle5_test_global_alternating_scores_df = as.data.frame(
  do.call(cbind, cycle5_test_global_alternating_scores_list))

cycle5_train_global_alternating_scores_max = apply(cycle5_train_global_alternating_scores_df, 1, 
                                                   function(row) {
                                                     return(max(unlist(row)))})
cycle5_train_global_alternating_scores_mean = apply(cycle5_train_global_alternating_scores_df, 1, 
                                                    function(row) {
                                                      return(mean(unlist(row)))})
cycle5_train_global_alternating_scores_min = apply(cycle5_train_global_alternating_scores_df, 1, 
                                                   function(row) {
                                                     return(min(unlist(row)))})

cycle5_train_global_alternating_scores_final = cbind(cycle5_train_global_alternating_scores_max,
                                                     cycle5_train_global_alternating_scores_mean,
                                                     cycle5_train_global_alternating_scores_min)
colnames(cycle5_train_global_alternating_scores_final) = paste0("global_alternating_scores_",
                                                                c("max", "mean", "min"))

saveRDS(cycle5_train_global_alternating_scores_final, "data/Created/tiling_global_alternating_scores_AorT_train.rds")


cycle5_test_global_alternating_scores_max = apply(cycle5_test_global_alternating_scores_df, 1, 
                                                  function(row) {
                                                    return(max(unlist(row)))})
cycle5_test_global_alternating_scores_mean = apply(cycle5_test_global_alternating_scores_df, 1, 
                                                   function(row) {
                                                     return(mean(unlist(row)))})
cycle5_test_global_alternating_scores_min = apply(cycle5_test_global_alternating_scores_df, 1, 
                                                  function(row) {
                                                    return(min(unlist(row)))})

cycle5_test_global_alternating_scores_final = cbind(cycle5_test_global_alternating_scores_max,
                                                    cycle5_test_global_alternating_scores_mean,
                                                    cycle5_test_global_alternating_scores_min)
colnames(cycle5_test_global_alternating_scores_final) = paste0("global_alternating_scores_",
                                                               c("max", "mean", "min"))

saveRDS(cycle5_test_global_alternating_scores_final, "data/Created/tiling_global_alternating_scores_AorT_test.rds")







# Random Library
cycle3_train = cycle3[random_train_indices,] %>%
  select(Sequence)
cycle3_test = cycle3[-random_train_indices,] %>%
  select(Sequence)
cycle3_all = rbind(cycle3_train, cycle3_test)

cycle3_pooling_by_5_scores_offset0 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=0)})
cycle3_pooling_by_5_scores_offset1 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=1)})
cycle3_pooling_by_5_scores_offset2 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=2)})
cycle3_pooling_by_5_scores_offset3 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=3)})
cycle3_pooling_by_5_scores_offset4 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=4)})

cycle3_all_pooling_by_5_scores_list = list(cycle3_pooling_by_5_scores_offset0,
                                             cycle3_pooling_by_5_scores_offset1,
                                             cycle3_pooling_by_5_scores_offset2,
                                             cycle3_pooling_by_5_scores_offset3,
                                             cycle3_pooling_by_5_scores_offset4)

cycle3_all_global_alternating_scores_list = lapply(cycle3_all_pooling_by_5_scores_list,
                                                     function(scores_list) {
                                                       lapply(scores_list, find_alternating_score)
                                                     })
cycle3_all_global_alternating_scores_df = as.data.frame(
  do.call(cbind, cycle3_all_global_alternating_scores_list))

cycle3_all_global_alternating_scores_max = apply(cycle3_all_global_alternating_scores_df, 1, 
                                                   function(row) {
                                                     return(max(unlist(row)))})
cycle3_all_global_alternating_scores_mean = apply(cycle3_all_global_alternating_scores_df, 1, 
                                                    function(row) {
                                                      return(mean(unlist(row)))})
cycle3_all_global_alternating_scores_min = apply(cycle3_all_global_alternating_scores_df, 1, 
                                                   function(row) {
                                                     return(min(unlist(row)))})

cycle3_all_global_alternating_scores_final = cbind(cycle3_all_global_alternating_scores_max,
                                                     cycle3_all_global_alternating_scores_mean,
                                                     cycle3_all_global_alternating_scores_min)
colnames(cycle3_all_global_alternating_scores_final) = paste0("global_alternating_scores_",
                                                                c("max", "mean", "min"))

saveRDS(cycle3_all_global_alternating_scores_final, "data/Created/random_all_global_alternating_scores_AorT.rds")








# A/T/AA/AT/TA/TT
cycle5_pooling_by_5_AorT_monodi_scores_offset0 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=0, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_pooling_by_5_AorT_monodi_scores_offset1 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=1, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_pooling_by_5_AorT_monodi_scores_offset2 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=2, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_pooling_by_5_AorT_monodi_scores_offset3 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=3, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_pooling_by_5_AorT_monodi_scores_offset4 = lapply(cycle5_train$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=4, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})

cycle5_train_pooling_by_5_AorT_monodi_scores_list = list(cycle5_pooling_by_5_AorT_monodi_scores_offset0,
                                             cycle5_pooling_by_5_AorT_monodi_scores_offset1,
                                             cycle5_pooling_by_5_AorT_monodi_scores_offset2,
                                             cycle5_pooling_by_5_AorT_monodi_scores_offset3,
                                             cycle5_pooling_by_5_AorT_monodi_scores_offset4)

cycle5_test_pooling_by_5_AorT_monodi_scores_offset0 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=0, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_test_pooling_by_5_AorT_monodi_scores_offset1 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=1, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_test_pooling_by_5_AorT_monodi_scores_offset2 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=2, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_test_pooling_by_5_AorT_monodi_scores_offset3 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=3, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle5_test_pooling_by_5_AorT_monodi_scores_offset4 = lapply(cycle5_test$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=4, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})

cycle5_test_pooling_by_5_AorT_monodi_scores_list = list(cycle5_test_pooling_by_5_AorT_monodi_scores_offset0,
                                            cycle5_test_pooling_by_5_AorT_monodi_scores_offset1,
                                            cycle5_test_pooling_by_5_AorT_monodi_scores_offset2,
                                            cycle5_test_pooling_by_5_AorT_monodi_scores_offset3,
                                            cycle5_test_pooling_by_5_AorT_monodi_scores_offset4)

cycle5_train_global_alternating_scores_AorT_monodi_list = lapply(cycle5_train_pooling_by_5_AorT_monodi_scores_list,
                                                     function(scores_list) {
                                                       lapply(scores_list, find_alternating_score)
                                                     })
cycle5_train_global_alternating_scores_AorT_monodi_df = as.data.frame(
  do.call(cbind, cycle5_train_global_alternating_scores_AorT_monodi_list))

cycle5_test_global_alternating_scores_AorT_monodi_list = lapply(cycle5_test_pooling_by_5_AorT_monodi_scores_list,
                                                    function(scores_list) {
                                                      lapply(scores_list, find_alternating_score)
                                                    })
cycle5_test_global_alternating_scores_AorT_monodi_df = as.data.frame(
  do.call(cbind, cycle5_test_global_alternating_scores_AorT_monodi_list))

cycle5_train_global_alternating_scores_AorT_monodi_max = apply(cycle5_train_global_alternating_scores_AorT_monodi_df, 1, 
                                                   function(row) {
                                                     return(max(unlist(row)))})
cycle5_train_global_alternating_scores_AorT_monodi_mean = apply(cycle5_train_global_alternating_scores_AorT_monodi_df, 1, 
                                                    function(row) {
                                                      return(mean(unlist(row)))})
cycle5_train_global_alternating_scores_AorT_monodi_min = apply(cycle5_train_global_alternating_scores_AorT_monodi_df, 1, 
                                                   function(row) {
                                                     return(min(unlist(row)))})

cycle5_train_global_alternating_scores_AorT_monodi_final = cbind(cycle5_train_global_alternating_scores_AorT_monodi_max,
                                                     cycle5_train_global_alternating_scores_AorT_monodi_mean,
                                                     cycle5_train_global_alternating_scores_AorT_monodi_min)
colnames(cycle5_train_global_alternating_scores_AorT_monodi_final) = paste0("global_alternating_scores_AorT_monodi_",
                                                                c("max", "mean", "min"))

saveRDS(cycle5_train_global_alternating_scores_AorT_monodi_final, "data/Created/tiling_global_alternating_scores_AorT_monodi_train.rds")


cycle5_test_global_alternating_scores_AorT_monodi_max = apply(cycle5_test_global_alternating_scores_AorT_monodi_df, 1, 
                                                  function(row) {
                                                    return(max(unlist(row)))})
cycle5_test_global_alternating_scores_AorT_monodi_mean = apply(cycle5_test_global_alternating_scores_AorT_monodi_df, 1, 
                                                   function(row) {
                                                     return(mean(unlist(row)))})
cycle5_test_global_alternating_scores_AorT_monodi_min = apply(cycle5_test_global_alternating_scores_AorT_monodi_df, 1, 
                                                  function(row) {
                                                    return(min(unlist(row)))})

cycle5_test_global_alternating_scores_AorT_monodi_final = cbind(cycle5_test_global_alternating_scores_AorT_monodi_max,
                                                    cycle5_test_global_alternating_scores_AorT_monodi_mean,
                                                    cycle5_test_global_alternating_scores_AorT_monodi_min)
colnames(cycle5_test_global_alternating_scores_AorT_monodi_final) = paste0("global_alternating_scores_AorT_monodi_",
                                                               c("max", "mean", "min"))

saveRDS(cycle5_test_global_alternating_scores_AorT_monodi_final, "data/Created/tiling_global_alternating_scores_AorT_monodi_test.rds")


# Random Library
cycle3_pooling_by_5_scores_AorT_monodi_offset0 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=0, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle3_pooling_by_5_scores_AorT_monodi_offset1 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=1, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle3_pooling_by_5_scores_AorT_monodi_offset2 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=2, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle3_pooling_by_5_scores_AorT_monodi_offset3 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=3, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})
cycle3_pooling_by_5_scores_AorT_monodi_offset4 = lapply(cycle3_all$Sequence, function(seq) {
  find_pooling_by_5_scores(seq, offset=4, substrs_of_interest = c("A","T","AA","AT","TA","TT"))})

cycle3_all_pooling_by_5_scores_AorT_monodi_list = list(cycle3_pooling_by_5_scores_AorT_monodi_offset0,
                                           cycle3_pooling_by_5_scores_AorT_monodi_offset1,
                                           cycle3_pooling_by_5_scores_AorT_monodi_offset2,
                                           cycle3_pooling_by_5_scores_AorT_monodi_offset3,
                                           cycle3_pooling_by_5_scores_AorT_monodi_offset4)

cycle3_all_global_alternating_scores_AorT_monodi_list = lapply(cycle3_all_pooling_by_5_scores_AorT_monodi_list,
                                                   function(scores_list) {
                                                     lapply(scores_list, find_alternating_score)
                                                   })
cycle3_all_global_alternating_scores_AorT_monodi_df = as.data.frame(
  do.call(cbind, cycle3_all_global_alternating_scores_AorT_monodi_list))

cycle3_all_global_alternating_scores_AorT_monodi_max = apply(cycle3_all_global_alternating_scores_AorT_monodi_df, 1, 
                                                 function(row) {
                                                   return(max(unlist(row)))})
cycle3_all_global_alternating_scores_AorT_monodi_mean = apply(cycle3_all_global_alternating_scores_AorT_monodi_df, 1, 
                                                  function(row) {
                                                    return(mean(unlist(row)))})
cycle3_all_global_alternating_scores_AorT_monodi_min = apply(cycle3_all_global_alternating_scores_AorT_monodi_df, 1, 
                                                 function(row) {
                                                   return(min(unlist(row)))})

cycle3_all_global_alternating_scores_AorT_monodi_final = cbind(cycle3_all_global_alternating_scores_AorT_monodi_max,
                                                   cycle3_all_global_alternating_scores_AorT_monodi_mean,
                                                   cycle3_all_global_alternating_scores_AorT_monodi_min)
colnames(cycle3_all_global_alternating_scores_AorT_monodi_final) = paste0("global_alternating_scores_AorT_monodi_",
                                                              c("max", "mean", "min"))

saveRDS(cycle3_all_global_alternating_scores_AorT_monodi_final, "data/Created/random_all_global_alternating_scores_AorT_monodi.rds")
