la.model_all_dist_freq_di_AorTdi_bid = readRDS("model/lasso-all_dist_freq_di_AorTdi_bid.rds")

# la.all_dist_freq_di_AorTdi_bid_best_lambda = la.model_all_dist_freq_di_AorTdi_bid$lambda.min
la.all_dist_freq_di_AorTdi_bid_best_lambda = la.model_all_dist_freq_di_AorTdi_bid$lambda.1se

la.all_dist_freq_di_AorTdi_bid_coefs = coef(la.model_all_dist_freq_di_AorTdi_bid, 
                                            s=la.all_dist_freq_di_AorTdi_bid_best_lambda)

#AA/AT/TA/TT Group
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and AA/AT/TA/TT Group") + 
  ylim(ylims)


# AC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and AA/AT/TA/TT Group") +
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TT and AA/AT/TA/TT Group") + 
  ylim(ylims)


alpha_val = .3
ggplot() +
  geom_point(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AorT"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Betas for Dist Freq Between All Dinucleotides and AA/AT/TA/TT Group") + 
  ylim(ylims)



#CC/CG/GC/GG Group:
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TT and CC/CG/GC/GG Group") + 
  ylim(ylims)


ggplot() +
  geom_point(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid_coefs[str_which(la.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CorG"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Betas for Dist Freq Between All Dinucleotides and CC/CG/GC/GG Group") + 
  ylim(ylims)





# Inputs: Position Specific Dinucleotides + Counts of Trinucleotides + 
# Fourier Features Based on all Dinucleotides + 
# Modulus of Difference of Fourier Features Based on all Dinucleotides with h=6 (period 9.8) Calculated using fft + 
# All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides

la.model_all_dist_freq_di_AorTdi_bid15 = readRDS("model/lasso-all_dist_freq_di_AorTdi_bid15.rds")

# la.all_dist_freq_di_AorTdi_bid15_best_lambda = la.model_all_dist_freq_di_AorTdi_bid15$lambda.min
la.all_dist_freq_di_AorTdi_bid15_best_lambda = la.model_all_dist_freq_di_AorTdi_bid15$lambda.1se

la.all_dist_freq_di_AorTdi_bid15_coefs = coef(la.model_all_dist_freq_di_AorTdi_bid15, 
                                              s=la.all_dist_freq_di_AorTdi_bid15_best_lambda)

#AA/AT/TA/TT Group
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TT and AA/AT/TA/TT Group") + 
  ylim(ylims)


ggplot() +
  geom_point(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AorT"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Betas for Dist Freq Between All Dinucleotides and AA/AT/TA/TT Group") + 
  ylim(ylims)


#CC/CG/GC/GG Group:
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "AT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "GT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "TT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TT and CC/CG/GC/GG Group") + 
  ylim(ylims)


ggplot() +
  geom_point(data=data.frame(Model.Coefficient = la.all_dist_freq_di_AorTdi_bid15_coefs[str_which(la.all_dist_freq_di_AorTdi_bid15_coefs@Dimnames[[1]], "CorG"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Betas for Dist Freq Between All Dinucleotides and CC/CG/GC/GG Group") + 
  ylim(ylims)





# Di Group lasso:

gla.model_all_dist_freq_di_AorTdi_bid = readRDS("model/group-lasso-all_dist_freq_di_AorTdi_bid.rds")

# gla.all_dist_freq_di_AorTdi_bid_best_lambda = gla.model_all_dist_freq_di_AorTdi_bid$lambda.min
gla.all_dist_freq_di_AorTdi_bid_best_lambda = gla.model_all_dist_freq_di_AorTdi_bid$lambda.1se

gla.all_dist_freq_di_AorTdi_bid_coefs = coef(gla.model_all_dist_freq_di_AorTdi_bid, 
                                             s=gla.all_dist_freq_di_AorTdi_bid_best_lambda)

#AA/AT/TA/TT Group
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and AA/AT/TA/TT Group") + 
  ylim(ylims)


# AC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and AA/AT/TA/TT Group") +
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TT and AA/AT/TA/TT Group") + 
  ylim(ylims)


alpha_val = .3
ggplot() +
  geom_point(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AorT"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Di group lasso Dist Freq Between All Dinucleotides and AA/AT/TA/TT Group") + 
  ylim(ylims)



#CC/CG/GC/GG Group:
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TT and CC/CG/GC/GG Group") + 
  ylim(ylims)


ggplot() +
  geom_point(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CorG"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Di group lasso Dist Freq Between All Dinucleotides and CC/CG/GC/GG Group") + 
  ylim(ylims)








# Dist Group lasso:

gla.model_all_dist_freq_di_AorTdi_bid = readRDS("model/group-lasso-all_dist_freq_di_AorTdi_bid.rds")

# gla.all_dist_freq_di_AorTdi_bid_best_lambda = gla.model_all_dist_freq_di_AorTdi_bid$lambda.min
gla.all_dist_freq_di_AorTdi_bid_best_lambda = gla.model_all_dist_freq_di_AorTdi_bid$lambda.1se

gla.all_dist_freq_di_AorTdi_bid_coefs = coef(gla.model_all_dist_freq_di_AorTdi_bid, 
                                             s=gla.all_dist_freq_di_AorTdi_bid_best_lambda)

#AA/AT/TA/TT Group
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and AA/AT/TA/TT Group") + 
  ylim(ylims)


# AC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and AA/AT/TA/TT Group") +
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TA_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TC_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TG_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TG and AA/AT/TA/TT Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TT_*_AorT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TT and AA/AT/TA/TT Group") + 
  ylim(ylims)


alpha_val = .3
ggplot() +
  geom_point(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AorT"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Di group lasso Dist Freq Between All Dinucleotides and AA/AT/TA/TT Group") + 
  ylim(ylims)



#CC/CG/GC/GG Group:
# AA:
ylims = c(-.05, .05)
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "AT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "GT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TA_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TC_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TG_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TG and CC/CG/GC/GG Group") + 
  ylim(ylims)

# TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "TT_*_CorG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TT and CC/CG/GC/GG Group") + 
  ylim(ylims)


ggplot() +
  geom_point(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_AorTdi_bid_coefs[str_which(gla.all_dist_freq_di_AorTdi_bid_coefs@Dimnames[[1]], "CorG"),],
                             Distance = rep(2:48, each=16)) %>% 
               mutate(Included.In.Model = (Model.Coefficient != 0)),
             aes(x=Distance,y=Model.Coefficient,color=Included.In.Model), alpha=alpha_val) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) + 
  ggtitle("Di group lasso Dist Freq Between All Dinucleotides and CC/CG/GC/GG Group") + 
  ylim(ylims)





