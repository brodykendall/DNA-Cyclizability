library(FactoMineR)

temp = FAMD(dat[,model_cols_all_int_inc_ratio.1])

temp.2 = MFA(dat[,model_cols_all_int_inc_ratio.1], group = c(50, 1),
             type = c("n", "s"))


library(multiway)


