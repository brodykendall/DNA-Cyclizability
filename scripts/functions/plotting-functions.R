construct_periodicity_plots_AT = function(data, var, ylims_quartiles=c(0.575, 0.65), 
                                          ylims_1000s=c(0.57, 0.77), include_title=TRUE,
                                          include_labels=TRUE, include_legend=TRUE) {
  
  var_string = deparse(substitute(var))
  data_string = deparse(substitute(data))
  var_of_interest = data %>% select({{var}})
  
  # Find cutoffs for first-fourth quartiles + custom
  cutoffs_0_25 = quantile(var_of_interest[,1], c(0, 0.25))
  cutoffs_25_50 = quantile(var_of_interest[,1], c(0.25, 0.5))
  cutoffs_50_75 = quantile(var_of_interest[,1], c(0.5, 0.75))
  cutoffs_75_100 = quantile(var_of_interest[,1], c(0.75, 1))
  cutoffs_custom1 = c(sort(var_of_interest[,1], TRUE)[1000], cutoffs_75_100[2])
  cutoffs_custom2 = c(cutoffs_0_25[1], sort(var_of_interest[,1], FALSE)[1000])
  
  # Divide each library into quartiles + custom
  var_q1 = data %>%
    filter({{var}} >= cutoffs_0_25[1] & {{var}} <= cutoffs_0_25[2]) %>%
    select(all_of(ps1), {{var}})
  var_q2 = data %>%
    filter({{var}} >= cutoffs_25_50[1] & {{var}} <= cutoffs_25_50[2]) %>%
    select(all_of(ps1), {{var}})
  var_q3 = data %>%
    filter({{var}} >= cutoffs_50_75[1] & {{var}} <= cutoffs_50_75[2]) %>%
    select(all_of(ps1), {{var}})
  var_q4 = data %>%
    filter({{var}} >= cutoffs_75_100[1] & {{var}} <= cutoffs_75_100[2]) %>%
    select(all_of(ps1), {{var}})
  var_custom1 = data %>%
    filter({{var}} >= cutoffs_custom1[1] & {{var}} <= cutoffs_custom1[2]) %>%
    select(all_of(ps1), {{var}})
  var_custom2 = data %>%
    filter({{var}} >= cutoffs_custom2[1] & {{var}} <= cutoffs_custom2[2]) %>%
    select(all_of(ps1), {{var}})
  
  # Find the relative frequencies of A and T at each position (1-50) for each quartile + custom
  A_T_var_q1 = apply(var_q1 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q1))
  })
  A_T_var_q2 = apply(var_q2 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q2))
  })
  A_T_var_q3 = apply(var_q3 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q3))
  })
  A_T_var_q4 = apply(var_q4 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q4))
  })
  A_T_var_custom1 = apply(var_custom1 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_custom1))
  })
  A_T_var_custom2 = apply(var_custom2 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_custom2))
  })
  
  # Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
  A_T_var_q = data.frame(q1 = A_T_var_q1,
                         q2 = A_T_var_q2,
                         q3 = A_T_var_q3,
                         q4 = A_T_var_q4, 
                         custom1 = A_T_var_custom1,
                         custom2 = A_T_var_custom2)
  
  # Plot the relative frequencies of A and T by position for each quartile (q1 v q4, 
  # q2 v q3, and custom)
  AT_q1vq4_var = ggplot(data = A_T_var_q, aes(x = 1:50)) +
    geom_line(aes(y=q1, color="0-25%")) +
    geom_line(aes(y=q4, color="75-100%")) +
    scale_color_manual(breaks = c("0-25%", "75-100%"),
                       values = c("blue", "red"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average A/T Content") +
    ylim(ylims_quartiles)
  
  AT_q2vq3_var = ggplot(data = A_T_var_q, aes(x = 1:50)) +
    geom_line(aes(y=q2, color="25-50%")) +
    geom_line(aes(y=q3, color="50-75%")) +
    scale_color_manual(breaks = c("25-50%", "50-75%"),
                       values = c("green", "orange"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average A/T Content") +
    ylim(ylims_quartiles)
  
  AT_custom1_var = ggplot(data = A_T_var_q, aes(x = 1:50)) +
    geom_line(aes(y=custom1, color="Top 1000")) +
    geom_line(aes(y=custom2, color="Bottom 1000")) +
    scale_color_manual(breaks = c("Top 1000", "Bottom 1000"),
                       values = c("firebrick1", "dodgerblue"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average A/T Content") +
    ylim(ylims_1000s)
  
  if(include_title) {
    AT_q1vq4_var = AT_q1vq4_var + labs(title = paste0("Average A/T Content by Position in ", data_string, ", ", var_string))
    AT_q2vq3_var = AT_q2vq3_var + labs(title = paste0("Average A/T Content by Position in ", data_string, ", ", var_string))
    AT_custom1_var = AT_custom1_var + labs(title = paste0("Average A/T Content by Position in ", data_string, ", ", var_string))
  }
  if(!include_labels) {
    AT_q1vq4_var = AT_q1vq4_var + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())
    AT_q2vq3_var = AT_q2vq3_var + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())
    AT_custom1_var = AT_custom1_var + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())
  }
  if(!include_legend) {
    AT_q1vq4_var = AT_q1vq4_var +
      theme(legend.position = "none")
    AT_q2vq3_var = AT_q2vq3_var +
      theme(legend.position = "none")
    AT_custom1_var = AT_custom1_var +
      theme(legend.position = "none")
  }
  
  plot_list = list(AT_q1vq4_var, AT_q2vq3_var, AT_custom1_var)
  
  return(plot_list)
}



construct_periodicity_plots_AA_TT_AT_TA = function(data, var, ylims_quartiles=c(0.35, 0.425), 
                                                   ylims_1000s=c(0.325, 0.65),
                                                   include_title=TRUE, include_labels=TRUE, 
                                                   include_legend=TRUE) {
  
  var_string = deparse(substitute(var))
  data_string = deparse(substitute(data))
  var_of_interest = data %>% select({{var}})
  
  # Find cutoffs for first-fourth quartiles + custom
  cutoffs_0_25 = quantile(var_of_interest[,1], c(0, 0.25))
  cutoffs_25_50 = quantile(var_of_interest[,1], c(0.25, 0.5))
  cutoffs_50_75 = quantile(var_of_interest[,1], c(0.5, 0.75))
  cutoffs_75_100 = quantile(var_of_interest[,1], c(0.75, 1))
  cutoffs_custom1 = c(sort(var_of_interest[,1], TRUE)[1000], cutoffs_75_100[2])
  cutoffs_custom2 = c(cutoffs_0_25[1], sort(var_of_interest[,1], FALSE)[1000])
  
  # Divide each library into quartiles + custom
  var_q1 = data %>%
    filter({{var}} >= cutoffs_0_25[1] & {{var}} <= cutoffs_0_25[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q2 = data %>%
    filter({{var}} >= cutoffs_25_50[1] & {{var}} <= cutoffs_25_50[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q3 = data %>%
    filter({{var}} >= cutoffs_50_75[1] & {{var}} <= cutoffs_50_75[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q4 = data %>%
    filter({{var}} >= cutoffs_75_100[1] & {{var}} <= cutoffs_75_100[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_custom1 = data %>%
    filter({{var}} >= cutoffs_custom1[1] & {{var}} <= cutoffs_custom1[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_custom2 = data %>%
    filter({{var}} >= cutoffs_custom2[1] & {{var}} <= cutoffs_custom2[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  
  # Find the relative frequencies of AA, TT, AT and TA at each position (1-49) for each quartile + custom
  AA_TT_AT_TA_var_q1 = apply(var_q1 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q1))
    # return((TA_freq)/nrow(var_q1))
  })
  AA_TT_AT_TA_var_q2 = apply(var_q2 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q2))
    # return((TA_freq)/nrow(var_q2))
  })
  AA_TT_AT_TA_var_q3 = apply(var_q3 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q3))
    # return((TA_freq)/nrow(var_q3))
  })
  AA_TT_AT_TA_var_q4 = apply(var_q4 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q4))
    # return((TA_freq)/nrow(var_q4))
  })
  AA_TT_AT_TA_var_custom1 = apply(var_custom1 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_custom1))
    # return((TA_freq)/nrow(var_custom1))
  })
  AA_TT_AT_TA_var_custom2 = apply(var_custom2 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_custom2))
    # return((TA_freq)/nrow(var_custom2))
  })
  
  # print("Note: Only TA_freq")
  
  # Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for all four quartiles + custom
  AA_TT_AT_TA_var_q = data.frame(q1 = AA_TT_AT_TA_var_q1,
                                 q2 = AA_TT_AT_TA_var_q2,
                                 q3 = AA_TT_AT_TA_var_q3,
                                 q4 = AA_TT_AT_TA_var_q4, 
                                 custom1 = AA_TT_AT_TA_var_custom1,
                                 custom2 = AA_TT_AT_TA_var_custom2)
  
  # Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4, 
  # q2 v q3, and custom)
  AATT_q1vq4_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
    geom_line(aes(y=q1, color="0-25%")) +
    geom_line(aes(y=q4, color="75-100%")) +
    scale_color_manual(breaks = c("0-25%", "75-100%"),
                       values = c("blue", "red"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average AA/TT/AT/TA Content") +
    ylim(ylims_quartiles)
  
  # AATT_q1vq4_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
  #   geom_line(aes(y=q1, color="0-25%")) + 
  #   geom_line(aes(y=q4, color="75-100%"))
  
  AATT_q2vq3_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
    geom_line(aes(y=q2, color="25-50%")) +
    geom_line(aes(y=q3, color="50-75%")) +
    scale_color_manual(breaks = c("25-50%", "50-75%"),
                       values = c("green", "orange"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average AA/TT/AT/TA Content") + 
    ylim(ylims_quartiles)
  
  AATT_custom1_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
    geom_line(aes(y=custom1, color="Top 1000")) +
    geom_line(aes(y=custom2, color="Bottom 1000")) +
    scale_color_manual(breaks = c("Top 1000", "Bottom 1000"),
                       values = c("firebrick1", "dodgerblue"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average AA/TT/AT/TA Content") +
    ylim(ylims_1000s)
  
  if(include_title) {
    AATT_q1vq4_var = AATT_q1vq4_var +
      labs(title = paste0("Average AA/TT/AT/TA Content by Position in ", data_string, ", ", var_string))
    AATT_q2vq3_var = AATT_q2vq3_var + 
      labs(title = paste0("Average AA/TT/AT/TA Content by Position in ", data_string, ", ", var_string))
    AATT_custom1_var = AATT_custom1_var + 
      labs(title = paste0("Average AA/TT/AT/TA Content by Position in ", data_string, ", ", var_string))
  }
  if(!include_labels) {
    AATT_q1vq4_var = AATT_q1vq4_var +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    AATT_q2vq3_var = AATT_q2vq3_var +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    AATT_custom1_var = AATT_custom1_var +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
  }
  if(!include_legend) {
    AATT_q1vq4_var = AATT_q1vq4_var +
      theme(legend.position = "none")
    AATT_q2vq3_var = AATT_q2vq3_var +
      theme(legend.position = "none")
    AATT_custom1_var = AATT_custom1_var +
      theme(legend.position = "none")
  }
  
  plot_list = list(AATT_q1vq4_var, AATT_q2vq3_var, AATT_custom1_var)
  
  return(plot_list)
}
