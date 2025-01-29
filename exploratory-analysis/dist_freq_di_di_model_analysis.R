la.model_all_dist_freq_di = readRDS("model/lasso-all_dist_freq_di.rds")

# la.all_dist_freq_di_AorTdi_bid_best_lambda = la.model_all_dist_freq_di_AorTdi_bid$lambda.min
la.all_dist_freq_di_best_lambda = la.model_all_dist_freq_di$lambda.1se

la.all_dist_freq_di_coefs = coef(la.model_all_dist_freq_di, 
                                 s=la.all_dist_freq_di_best_lambda)

ylims = c(-.055, .055)
# AA_AA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and AA") + 
  ylim(ylims)

# AA_AC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and AC") + 
  ylim(ylims)

# AA_AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and AG") + 
  ylim(ylims)

# AA_AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and AT") + 
  ylim(ylims)

# AA_CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and CA") + 
  ylim(ylims)

# AA_CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and CC") + 
  ylim(ylims)

# AA_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and CG") + 
  ylim(ylims)

# AA_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and CT") + 
  ylim(ylims)

# AA_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and GA") + 
  ylim(ylims)

# AA_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and GC") + 
  ylim(ylims)

# AA_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and GG") + 
  ylim(ylims)

# AA_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and GT") + 
  ylim(ylims)

# AA_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and TA") + 
  ylim(ylims)

# AA_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and TC") + 
  ylim(ylims)

# AA_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and TG") + 
  ylim(ylims)

# AA_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AA and TT") + 
  ylim(ylims)


# AC_AC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_AC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and AC") + 
  ylim(ylims)

# AC_AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_AG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and AG") + 
  ylim(ylims)

# AC_AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and AT") + 
  ylim(ylims)

# AC_CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and CA") + 
  ylim(ylims)

# AC_CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and CC") + 
  ylim(ylims)

# AC_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and CG") + 
  ylim(ylims)

# AC_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and CT") + 
  ylim(ylims)

# AC_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and GA") + 
  ylim(ylims)

# AC_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and GC") + 
  ylim(ylims)

# AC_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and GG") + 
  ylim(ylims)

# AC_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and GT") + 
  ylim(ylims)

# AC_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and TA") + 
  ylim(ylims)

# AC_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and TC") + 
  ylim(ylims)

# AC_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and TG") + 
  ylim(ylims)

# AC_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AC and TT") + 
  ylim(ylims)


# AG_AG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_AG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and AG") + 
  ylim(ylims)

# AG_AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and AT") + 
  ylim(ylims)

# AG_CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and CA") + 
  ylim(ylims)

# AG_CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and CC") + 
  ylim(ylims)

# AG_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and CG") + 
  ylim(ylims)

# AG_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and CT") + 
  ylim(ylims)

# AG_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and GA") + 
  ylim(ylims)

# AG_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and GC") + 
  ylim(ylims)

# AG_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and GG") + 
  ylim(ylims)

# AG_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and GT") + 
  ylim(ylims)

# AG_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and TA") + 
  ylim(ylims)

# AG_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and TC") + 
  ylim(ylims)

# AG_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and TG") + 
  ylim(ylims)

# AG_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AG and TT") + 
  ylim(ylims)


# AT_AT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and AT") + 
  ylim(ylims)

# AT_CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and CA") + 
  ylim(ylims)

# AT_CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and CC") + 
  ylim(ylims)

# AT_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and CG") + 
  ylim(ylims)

# AT_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and CT") + 
  ylim(ylims)

# AT_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and GA") + 
  ylim(ylims)

# AT_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and GC") + 
  ylim(ylims)

# AT_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and GG") + 
  ylim(ylims)

# AT_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and GT") + 
  ylim(ylims)

# AT_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and TA") + 
  ylim(ylims)

# AT_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and TC") + 
  ylim(ylims)

# AT_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and TG") + 
  ylim(ylims)

# AT_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between AT and TT") + 
  ylim(ylims)


# CA_CA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and CA") + 
  ylim(ylims)

# CA_CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and CC") + 
  ylim(ylims)

# CA_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and CG") + 
  ylim(ylims)

# CA_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and CT") + 
  ylim(ylims)

# CA_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and GA") + 
  ylim(ylims)

# CA_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and GC") + 
  ylim(ylims)

# CA_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and GG") + 
  ylim(ylims)

# CA_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and GT") + 
  ylim(ylims)

# CA_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and TA") + 
  ylim(ylims)

# CA_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and TC") + 
  ylim(ylims)

# CA_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and TG") + 
  ylim(ylims)

# CA_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CA and TT") + 
  ylim(ylims)


# CC_CC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and CC") + 
  ylim(ylims)

# CC_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and CG") + 
  ylim(ylims)

# CC_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and CT") + 
  ylim(ylims)

# CC_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and GA") + 
  ylim(ylims)

# CC_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and GC") + 
  ylim(ylims)

# CC_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and GG") + 
  ylim(ylims)

# CC_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and GT") + 
  ylim(ylims)

# CC_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and TA") + 
  ylim(ylims)

# CC_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and TC") + 
  ylim(ylims)

# CC_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and TG") + 
  ylim(ylims)

# CC_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CC and TT") + 
  ylim(ylims)


# CG_CG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and CG") + 
  ylim(ylims)

# CG_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and CT") + 
  ylim(ylims)

# CG_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and GA") + 
  ylim(ylims)

# CG_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and GC") + 
  ylim(ylims)

# CG_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and GG") + 
  ylim(ylims)

# CG_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and GT") + 
  ylim(ylims)

# CG_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and TA") + 
  ylim(ylims)

# CG_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and TC") + 
  ylim(ylims)

# CG_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and TG") + 
  ylim(ylims)

# CG_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CG and TT") + 
  ylim(ylims)


# CT_CT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and CT") + 
  ylim(ylims)

# CT_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and GA") + 
  ylim(ylims)

# CT_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and GC") + 
  ylim(ylims)

# CT_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and GG") + 
  ylim(ylims)

# CT_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and GT") + 
  ylim(ylims)

# CT_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and TA") + 
  ylim(ylims)

# CT_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and TC") + 
  ylim(ylims)

# CT_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and TG") + 
  ylim(ylims)

# CT_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between CT and TT") + 
  ylim(ylims)


# GA_GA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and GA") + 
  ylim(ylims)

# GA_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and GC") + 
  ylim(ylims)

# GA_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and GG") + 
  ylim(ylims)

# GA_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and GT") + 
  ylim(ylims)

# GA_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and TA") + 
  ylim(ylims)

# GA_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and TC") + 
  ylim(ylims)

# GA_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and TG") + 
  ylim(ylims)

# GA_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GA and TT") + 
  ylim(ylims)


# GC_GC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and GC") + 
  ylim(ylims)

# GC_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and GG") + 
  ylim(ylims)

# GC_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and GT") + 
  ylim(ylims)

# GC_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and TA") + 
  ylim(ylims)

# GC_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and TC") + 
  ylim(ylims)

# GC_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and TG") + 
  ylim(ylims)

# GC_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GC and TT") + 
  ylim(ylims)


# GG_GG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GG_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and GG") + 
  ylim(ylims)

# GG_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GG_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and GT") + 
  ylim(ylims)

# GG_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and TA") + 
  ylim(ylims)

# GG_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and TC") + 
  ylim(ylims)

# GG_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and TG") + 
  ylim(ylims)

# GG_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GG and TT") + 
  ylim(ylims)


# GT_GT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GT_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and GT") + 
  ylim(ylims)

# GT_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and TA") + 
  ylim(ylims)

# GT_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and TC") + 
  ylim(ylims)

# GT_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and TG") + 
  ylim(ylims)

# GT_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between GT and TT") + 
  ylim(ylims)


# TA_TA:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and TA") + 
  ylim(ylims)

# TA_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and TC") + 
  ylim(ylims)

# TA_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and TG") + 
  ylim(ylims)

# TA_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TA and TT") + 
  ylim(ylims)


# TC_TC:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and TC") + 
  ylim(ylims)

# TC_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and TG") + 
  ylim(ylims)

# TC_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TC and TT") + 
  ylim(ylims)


# TG_TG:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TG and TG") + 
  ylim(ylims)

# TG_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TG and TT") + 
  ylim(ylims)


# TT_TT:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "TT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Betas for Distance Frequency Between TT and TT") + 
  ylim(ylims)



# Similar plots based on distances rather than pairs of dinucleotides:

dinucleotides_combinations_regexes = as.matrix(expand.grid(di1=dinucleotides, di2=dinucleotides) %>%
                                                 filter(as.character(di1) <= as.character(di2))) %>%
  apply(1, function(row) {
    return(paste(row[1], row[2], sep="_"))
  })

dinucleotides_combinations_df = expand.grid(di1=dinucleotides, di2=dinucleotides) %>%
  filter(as.character(di1) <= as.character(di2))

# dist2:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist2$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist2")

# dist3:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist3$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist3")

# dist4:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist4$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist4")

# dist5:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist5$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist5")

# dist6:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist6$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist6")

# dist7:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist7$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist7")

# dist8:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist8$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist8")

# dist9:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist9$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist9")

# dist10:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist10$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist10")

# dist11:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist11$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist11")

# dist12:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist12$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist12")

# dist13:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist13$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist13")

# dist14:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist14$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist14")

# dist15:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist15$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist15")

# dist16:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist16$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist16")

# dist17:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist17$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist17")

# dist18:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist18$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist18")

# dist19:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist19$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist19")

# dist20:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist20$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist20")

# dist21:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist21$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist21")

# dist22:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist22$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist22")

# dist23:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist23$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist23")

# dist24:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist24$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist24")

# dist25:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist25$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist25")

# dist26:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist26$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist26")

# dist27:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist27$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist27")

# dist28:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist28$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist28")

# dist29:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist29$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist29")

# dist30:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist30$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist30")

# dist31:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist31$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist31")

# dist32:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist32$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist32")

# dist33:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist33$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist33")

# dist34:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist34$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist34")

# dist35:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist35$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist35")

# dist36:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist36$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist36")

# dist37:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist37$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist37")

# dist38:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist38$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist38")

# dist39:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist39$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist39")

# dist40:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist40$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist40")

# dist41:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist41$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist41")

# dist42:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist42$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist42")

# dist43:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist43$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist43")

# dist44:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist44$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist44")

# dist45:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist45$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist45")

# dist46:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist46$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist46")

# dist47:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist47$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist47")

# dist48:
ggplot(data=data.frame(Model.Coefficient = la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], "dist48$"),],
                       Dinucleotide1 = dinucleotides_combinations_df$di1,
                       Dinucleotide2 = dinucleotides_combinations_df$di2) %>% 
         mutate(Included.In.Model = factor((Model.Coefficient != 0)))) +
  geom_tile(aes(x=Dinucleotide1, y=Dinucleotide2, group=Model.Coefficient, fill=Model.Coefficient, alpha=Included.In.Model)) +
  scale_alpha_manual(values=c("TRUE"=1 ,"FALSE"=0)) +
  scale_fill_gradient(limits = ylims) +
  ggtitle("Betas for Dinucleotide Pairs at dist48")

# Try constructing model(s) with subsets of the features based on:
# (1) how many in a certain group are nonzero
# (2) whether the sum of the absolute value of the coefficients is above a certain threshold
# (3) some combination of (1) and (2)

# (I feel like you can generally look at the plots above and determine which 
# should be included and which shouldn't, just need to find a good metric to capture this)



# (1)

sufficient_nonzero = map(dinucleotides_combinations_regexes, ~la.all_dist_freq_di_coefs[str_which(la.all_dist_freq_di_coefs@Dimnames[[1]], .x),]) %>%
  map(~sum(.x==0) < 24) %>%
  unlist()

names(sufficient_nonzero) = dinucleotides_combinations_regexes









# Di Group lasso:

gla.model_all_dist_freq_di = readRDS("model/group-lasso-all_dist_freq_di.rds")

# gla.all_dist_freq_di_AorTdi_bid_best_lambda = gla.model_all_dist_freq_di_AorTdi_bid$lambda.min
gla.all_dist_freq_di_best_lambda = gla.model_all_dist_freq_di$lambda.1se

gla.all_dist_freq_di_coefs = coef(gla.model_all_dist_freq_di, 
                                 s=gla.all_dist_freq_di_best_lambda)

# AA_AA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and AA") + 
  ylim(ylims)

# AA_AC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and AC") + 
  ylim(ylims)

# AA_AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and AG") + 
  ylim(ylims)

# AA_AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and AT") + 
  ylim(ylims)

# AA_CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and CA") + 
  ylim(ylims)

# AA_CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and CC") + 
  ylim(ylims)

# AA_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and CG") + 
  ylim(ylims)

# AA_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and CT") + 
  ylim(ylims)

# AA_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and GA") + 
  ylim(ylims)

# AA_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and GC") + 
  ylim(ylims)

# AA_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and GG") + 
  ylim(ylims)

# AA_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and GT") + 
  ylim(ylims)

# AA_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and TA") + 
  ylim(ylims)

# AA_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and TC") + 
  ylim(ylims)

# AA_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and TG") + 
  ylim(ylims)

# AA_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AA and TT") + 
  ylim(ylims)


# AC_AC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_AC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and AC") + 
  ylim(ylims)

# AC_AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_AG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and AG") + 
  ylim(ylims)

# AC_AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and AT") + 
  ylim(ylims)

# AC_CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and CA") + 
  ylim(ylims)

# AC_CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and CC") + 
  ylim(ylims)

# AC_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and CG") + 
  ylim(ylims)

# AC_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and CT") + 
  ylim(ylims)

# AC_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and GA") + 
  ylim(ylims)

# AC_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and GC") + 
  ylim(ylims)

# AC_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and GG") + 
  ylim(ylims)

# AC_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and GT") + 
  ylim(ylims)

# AC_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and TA") + 
  ylim(ylims)

# AC_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and TC") + 
  ylim(ylims)

# AC_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and TG") + 
  ylim(ylims)

# AC_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AC and TT") + 
  ylim(ylims)


# AG_AG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_AG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and AG") + 
  ylim(ylims)

# AG_AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and AT") + 
  ylim(ylims)

# AG_CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and CA") + 
  ylim(ylims)

# AG_CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and CC") + 
  ylim(ylims)

# AG_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and CG") + 
  ylim(ylims)

# AG_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and CT") + 
  ylim(ylims)

# AG_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and GA") + 
  ylim(ylims)

# AG_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and GC") + 
  ylim(ylims)

# AG_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and GG") + 
  ylim(ylims)

# AG_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and GT") + 
  ylim(ylims)

# AG_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and TA") + 
  ylim(ylims)

# AG_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and TC") + 
  ylim(ylims)

# AG_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and TG") + 
  ylim(ylims)

# AG_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AG and TT") + 
  ylim(ylims)


# AT_AT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_AT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and AT") + 
  ylim(ylims)

# AT_CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and CA") + 
  ylim(ylims)

# AT_CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and CC") + 
  ylim(ylims)

# AT_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and CG") + 
  ylim(ylims)

# AT_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and CT") + 
  ylim(ylims)

# AT_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and GA") + 
  ylim(ylims)

# AT_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and GC") + 
  ylim(ylims)

# AT_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and GG") + 
  ylim(ylims)

# AT_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and GT") + 
  ylim(ylims)

# AT_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and TA") + 
  ylim(ylims)

# AT_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and TC") + 
  ylim(ylims)

# AT_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and TG") + 
  ylim(ylims)

# AT_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "AT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between AT and TT") + 
  ylim(ylims)


# CA_CA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and CA") + 
  ylim(ylims)

# CA_CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and CC") + 
  ylim(ylims)

# CA_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and CG") + 
  ylim(ylims)

# CA_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and CT") + 
  ylim(ylims)

# CA_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and GA") + 
  ylim(ylims)

# CA_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and GC") + 
  ylim(ylims)

# CA_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and GG") + 
  ylim(ylims)

# CA_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and GT") + 
  ylim(ylims)

# CA_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and TA") + 
  ylim(ylims)

# CA_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and TC") + 
  ylim(ylims)

# CA_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and TG") + 
  ylim(ylims)

# CA_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CA and TT") + 
  ylim(ylims)


# CC_CC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_CC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and CC") + 
  ylim(ylims)

# CC_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and CG") + 
  ylim(ylims)

# CC_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and CT") + 
  ylim(ylims)

# CC_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and GA") + 
  ylim(ylims)

# CC_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and GC") + 
  ylim(ylims)

# CC_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and GG") + 
  ylim(ylims)

# CC_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and GT") + 
  ylim(ylims)

# CC_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and TA") + 
  ylim(ylims)

# CC_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and TC") + 
  ylim(ylims)

# CC_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and TG") + 
  ylim(ylims)

# CC_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CC and TT") + 
  ylim(ylims)


# CG_CG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_CG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and CG") + 
  ylim(ylims)

# CG_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and CT") + 
  ylim(ylims)

# CG_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and GA") + 
  ylim(ylims)

# CG_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and GC") + 
  ylim(ylims)

# CG_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and GG") + 
  ylim(ylims)

# CG_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and GT") + 
  ylim(ylims)

# CG_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and TA") + 
  ylim(ylims)

# CG_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and TC") + 
  ylim(ylims)

# CG_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and TG") + 
  ylim(ylims)

# CG_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CG and TT") + 
  ylim(ylims)


# CT_CT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_CT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and CT") + 
  ylim(ylims)

# CT_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and GA") + 
  ylim(ylims)

# CT_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and GC") + 
  ylim(ylims)

# CT_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and GG") + 
  ylim(ylims)

# CT_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and GT") + 
  ylim(ylims)

# CT_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and TA") + 
  ylim(ylims)

# CT_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and TC") + 
  ylim(ylims)

# CT_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and TG") + 
  ylim(ylims)

# CT_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "CT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between CT and TT") + 
  ylim(ylims)


# GA_GA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and GA") + 
  ylim(ylims)

# GA_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and GC") + 
  ylim(ylims)

# GA_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and GG") + 
  ylim(ylims)

# GA_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and GT") + 
  ylim(ylims)

# GA_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and TA") + 
  ylim(ylims)

# GA_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and TC") + 
  ylim(ylims)

# GA_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and TG") + 
  ylim(ylims)

# GA_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GA and TT") + 
  ylim(ylims)


# GC_GC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_GC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and GC") + 
  ylim(ylims)

# GC_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and GG") + 
  ylim(ylims)

# GC_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and GT") + 
  ylim(ylims)

# GC_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and TA") + 
  ylim(ylims)

# GC_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and TC") + 
  ylim(ylims)

# GC_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and TG") + 
  ylim(ylims)

# GC_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GC and TT") + 
  ylim(ylims)


# GG_GG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GG_GG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and GG") + 
  ylim(ylims)

# GG_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GG_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and GT") + 
  ylim(ylims)

# GG_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and TA") + 
  ylim(ylims)

# GG_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and TC") + 
  ylim(ylims)

# GG_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and TG") + 
  ylim(ylims)

# GG_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GG and TT") + 
  ylim(ylims)


# GT_GT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GT_GT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and GT") + 
  ylim(ylims)

# GT_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and TA") + 
  ylim(ylims)

# GT_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and TC") + 
  ylim(ylims)

# GT_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and TG") + 
  ylim(ylims)

# GT_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "GT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between GT and TT") + 
  ylim(ylims)


# TA_TA:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TA"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and TA") + 
  ylim(ylims)

# TA_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and TC") + 
  ylim(ylims)

# TA_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and TG") + 
  ylim(ylims)

# TA_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TA_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TA and TT") + 
  ylim(ylims)


# TC_TC:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TC_TC"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and TC") + 
  ylim(ylims)

# TC_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TC_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and TG") + 
  ylim(ylims)

# TC_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TC_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TC and TT") + 
  ylim(ylims)


# TG_TG:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TG_TG"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TG and TG") + 
  ylim(ylims)

# TG_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TG_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TG and TT") + 
  ylim(ylims)


# TT_TT:
ggplot(data=data.frame(Model.Coefficient = gla.all_dist_freq_di_coefs[str_which(gla.all_dist_freq_di_coefs@Dimnames[[1]], "TT_TT"),],
                       Distance = 2:48) %>% 
         mutate(Included.In.Model = (Model.Coefficient != 0))) +
  geom_point(aes(x=Distance,y=Model.Coefficient,color=Included.In.Model)) + 
  scale_color_manual(values=c("TRUE" = "black", "FALSE" = "red")) + 
  geom_smooth(aes(x=Distance, y=Model.Coefficient),span = .4, level=0, color="black", size=.7) +
  ggtitle("Di group lasso Distance Frequency Between TT and TT") + 
  ylim(ylims)




