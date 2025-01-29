# Group lasso
library(gglasso)

# or:
library(grplasso)

# or:
# This one should be best as we (should) have sparsity within groups and between groups
# https://arxiv.org/pdf/2208.02942.pdf
library(sparsegl)


# Inputs: All Distance Frequency Between AA/TT/AT/TA and CC/GG/CG/GC Dinucleotides Only

X_all_dist_freq_AorTdi_bid = dat_all_dist_freq_AorTdi_bid %>%
  select(-C0_new)

## Groupings based on pairs of Dinucleotide Groups (3 groups)
ggla.di_groups_all_dist_freq_AorTdi_bid = rep(1:3, each=47)

set.seed(50)
ggla.all_dist_freq_AorTdi_bid = cv.gglasso(as.matrix(X_all_dist_freq_AorTdi_bid), y, group=ggla.di_groups_all_dist_freq_AorTdi_bid)

cor(all_dist_freq_AorTdi_bid_lm$fitted.values, y)