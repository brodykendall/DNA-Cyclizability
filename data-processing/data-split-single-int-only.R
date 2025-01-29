library(tidyverse)
library(modelr)

set.seed(50)

dat <- readRDS("data/Created/processed_single_int_only.rds")

kfold_obj <- dat %>%
  crossv_kfold(10, id = "fold") %>%
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame))

saveRDS(kfold_obj, "data/Created/train-10-fold-single-int-only.rds")