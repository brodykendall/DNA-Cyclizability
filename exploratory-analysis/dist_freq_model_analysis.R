library(glmnet)

### Regression on Dist Freq Features for Specific Dinucleotides:

# AA:
dat_all_dist_freq_AA_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("AA")),
                                        dat %>% select(C0_new))

all_dist_freq_AA_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_AA_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AA_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between AA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AC:
dat_all_dist_freq_AC_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("AC")),
                                        dat %>% select(C0_new))

all_dist_freq_AC_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_AC_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AC_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between AC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AG:
dat_all_dist_freq_AG_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("AG")),
                                        dat %>% select(C0_new))

all_dist_freq_AG_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_AG_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AG_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between AG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AT:
dat_all_dist_freq_AT_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("AT")),
                                        dat %>% select(C0_new))

all_dist_freq_AT_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_AT_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AT_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between AT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CA:
dat_all_dist_freq_CA_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("CA")),
                                        dat %>% select(C0_new))

all_dist_freq_CA_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_CA_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_CA_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between CA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CC:
dat_all_dist_freq_CC_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("CC")),
                                        dat %>% select(C0_new))

all_dist_freq_CC_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_CC_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_CC_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between CC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CG:
dat_all_dist_freq_CG_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("CG")),
                                        dat %>% select(C0_new))

all_dist_freq_CG_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_CG_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_CG_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between CG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CT:
dat_all_dist_freq_CT_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("CT")),
                                        dat %>% select(C0_new))

all_dist_freq_CT_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_CT_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_CT_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between CT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GA:
dat_all_dist_freq_GA_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("GA")),
                                        dat %>% select(C0_new))

all_dist_freq_GA_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_GA_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_GA_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between GA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GC:
dat_all_dist_freq_GC_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("GC")),
                                        dat %>% select(C0_new))

all_dist_freq_GC_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_GC_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_GC_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between GC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GG:
dat_all_dist_freq_GG_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("GG")),
                                        dat %>% select(C0_new))

all_dist_freq_GG_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_GG_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_GG_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between GG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GT:
dat_all_dist_freq_GT_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("GT")),
                                        dat %>% select(C0_new))

all_dist_freq_GT_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_GT_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_GT_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between GT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TA:
dat_all_dist_freq_TA_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("TA")),
                                        dat %>% select(C0_new))

all_dist_freq_TA_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_TA_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_TA_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between TA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TC:
dat_all_dist_freq_TC_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("TC")),
                                        dat %>% select(C0_new))

all_dist_freq_TC_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_TC_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_TC_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between TC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TG:
dat_all_dist_freq_TG_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("TG")),
                                        dat %>% select(C0_new))

all_dist_freq_TG_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_TG_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_TG_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between TG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TT:
dat_all_dist_freq_TT_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% select(starts_with("TT")),
                                        dat %>% select(C0_new))

all_dist_freq_TT_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_TT_AorTdi_bid)


ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_TT_AorTdi_bid_lm)[-1],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Reduced Model Betas for Dist Freq Between TT and AA/AT/TA/TT Group") + 
  ylim(ylims)




### Regression on Dist Freq Features for AA, AT, TA, and TT:

dat_all_dist_freq_AAAT_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% 
                                            select(starts_with("AA"), starts_with("AT"),
                                                   starts_with("TA"), starts_with("TT")),
                                          dat %>% select(C0_new))

all_dist_freq_AAAT_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_AAAT_AorTdi_bid)


# AA:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AAAT_AorTdi_bid_lm)[2:48],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("AA Betas. Model includes AA, AT, TA, and TT") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AAAT_AorTdi_bid_lm)[49:95],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("AT Betas. Model includes AA, AT, TA, and TT") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AAAT_AorTdi_bid_lm)[96:142],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("TA Betas. Model includes AA, AT, TA, and TT") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_AAAT_AorTdi_bid_lm)[143:189],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("TT Betas. Model includes AA, AT, TA, and TT") + 
  ylim(ylims)


### Regression on Dist Freq Features for AA, AC, AG, and AT:

dat_all_dist_freq_A_AorTdi_bid = cbind(dist_freq_di_AorTdi_bid2 %>% 
                                         select(starts_with("AA"), starts_with("AC"),
                                                starts_with("AG"), starts_with("AT")),
                                       dat %>% select(C0_new))

all_dist_freq_A_AorTdi_bid_lm = lm(C0_new~., data=dat_all_dist_freq_A_AorTdi_bid)


# AA:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_A_AorTdi_bid_lm)[2:48],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("AA Betas. Model includes AA, AC, AG, and AT") + 
  ylim(ylims)

# AC:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_A_AorTdi_bid_lm)[49:95],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("AC Betas. Model includes AA, AC, AG, and AT") + 
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_A_AorTdi_bid_lm)[96:142],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("AG Betas. Model includes AA, AC, AG, and AT") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = coef(all_dist_freq_A_AorTdi_bid_lm)[143:189],
                       Distance = 2:48)) +
  geom_point(aes(x=Distance,y=Model.Coefficient), color="black") + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("AT Betas. Model includes AA, AC, AG, and AT") + 
  ylim(ylims)