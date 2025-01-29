normalize_column <- function(column, max_val, min_val) {
  scaled_column <- (column - min_val) / (max_val - min_val)
  scaled_column <- scaled_column * 2 - 1
  return(scaled_column)
}

tiling_max_val = max(tiling %>% select(C26, C29, C31, C0))
tiling_min_val = min(tiling %>% select(C26, C29, C31, C0))

tiling = tiling %>% mutate(
  C26_scaled = normalize_column(C26, tiling_max_val, tiling_min_val),
  C29_scaled = normalize_column(C29, tiling_max_val, tiling_min_val),
  C31_scaled = normalize_column(C31, tiling_max_val, tiling_min_val),
  C0_scaled = normalize_column(C0, tiling_max_val, tiling_min_val),
  amp_scaled = normalize_column(Amplitude, tiling_max_val, tiling_min_val))

chrv_max_val = max(chrv %>% select(C26, C29, C31, C0))
chrv_min_val = min(chrv %>% select(C26, C29, C31, C0))

chrv = chrv %>% mutate(
  C26_scaled = normalize_column(C26, chrv_max_val, chrv_min_val),
  C29_scaled = normalize_column(C29, chrv_max_val, chrv_min_val),
  C31_scaled = normalize_column(C31, chrv_max_val, chrv_min_val),
  C0_scaled = normalize_column(C0, chrv_max_val, chrv_min_val),
  amp_scaled = normalize_column(Amplitude, chrv_max_val, chrv_min_val))

plot_C_by_library_position = function(data, starting_point=1, num_sequences=10, include_C0=FALSE) {
  row_sequence = seq(from=starting_point, by=1, length.out=num_sequences)
  position_sequence = (row_sequence-1)*7
  # sin_position_sequence = seq((starting_point-1)*7, (starting_point-1)*7+7*(num_sequences-1))
  
  if (include_C0) { 
    alpha1=1
    alpha2=0.7
    
    ggplot(mapping=aes(x=position_sequence, y=data[row_sequence,])) +
      # geom_line(aes(y=data[row_sequence, "C26"], color="C26"), alpha=alpha1) +
      # geom_line(aes(y=data[row_sequence, "C29"], color="C29"), alpha=alpha1) + 
      # geom_line(aes(y=data[row_sequence, "C31"], color="C31"), alpha=alpha1) +
      geom_point(aes(y=data[row_sequence, "C26"], color="C26"), alpha=alpha1) +
      geom_point(aes(y=data[row_sequence, "C29"], color="C29"), alpha=alpha1) + 
      geom_point(aes(y=data[row_sequence, "C31"], color="C31"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence,"C0"], color="C0"), alpha=alpha2) +
      # geom_line(aes(y=apply(data[row_sequence,], 1,
      #                       function(row) {mean(as.numeric(row[c("C26", "C31")]))}),
      #               color="mean C26, C31")) +
      # geom_line(aes(x=sin_position_sequence,
      #               y=sin_amp_multiplier*sin(2*pi*sin_position_sequence/10.4), color="sin"), alpha=0.3) +
      # geom_line(aes(y=apply(data[row_sequence,], 1, 
      #                       function(row) {mean(as.numeric(row[c("C26", "C31")]))}) + 
      #                 sin_amp_multiplier*sin(2*pi*position_sequence/10.4), color="mean C26, C31 - sin")) +
      xlab("Position") +
      ylab("Value") +
      scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green", "C0"="orange"))
  }
  else {
    alpha1=1
    # alpha2=1
    
    ggplot(mapping=aes(x=position_sequence, y=data[row_sequence,])) +
      geom_line(aes(y=data[row_sequence, "C26"], color="C26"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence, "C29"], color="C29"), alpha=alpha1) + 
      geom_line(aes(y=data[row_sequence, "C31"], color="C31"), alpha=alpha1) +
      # geom_line(aes(y=data[row_sequence,"C0"], color="C0"), alpha=alpha2) +
      # geom_line(aes(y=apply(data[row_sequence,], 1,
      #                       function(row) {mean(as.numeric(row[c("C26", "C31")]))}),
      #               color="mean C26, C31")) +
      # geom_line(aes(x=sin_position_sequence,
      #               y=sin_amp_multiplier*sin(2*pi*sin_position_sequence/10.4), color="sin"), alpha=0.3) +
      # geom_line(aes(y=apply(data[row_sequence,], 1, 
      #                       function(row) {mean(as.numeric(row[c("C26", "C31")]))}) + 
      #                 sin_amp_multiplier*sin(2*pi*position_sequence/10.4), color="mean C26, C31 - sin")) +
      xlab("Position") +
      ylab("Value") +
      scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green"))
  }
}

plot_rankings_C_by_library_position = function(data, center_point=2, num_sequences=3, aa=c(1,1/0.82,1/0.7),
                                               C0_var = C0, phase_var = Phase, amp_var = Amplitude) {
  row_sequence = seq(from=(center_point-(num_sequences+1)/2), by=1, length.out=num_sequences)
  position_sequence = (row_sequence-1)*7
  
  data = data[row_sequence,]
  
  
  
  # seq_of_interest = (num_sequences+1)/2
  # C0_of_interest = data$C0_scaled[seq_of_interest]
  # phase_of_interest = data$Phase[seq_of_interest]
  # amp_of_interest = data$amp_scaled[seq_of_interest]
  seq_of_interest = (num_sequences+1)/2
  # C0_of_interest = data$C0[seq_of_interest]
  C0_of_interest = (data %>% select({{C0_var}}))[seq_of_interest,][[1]]
  # phase_of_interest = data$Phase[seq_of_interest]
  phase_of_interest = (data %>% select({{phase_var}}))[seq_of_interest,][[1]]
  # amp_of_interest = data$Amplitude[seq_of_interest]
  amp_of_interest = (data %>% select({{amp_var}}))[seq_of_interest,][[1]]
  
  # print((data %>% select(C26_scaled, C29_scaled, C31_scaled, C0_scaled, amp_scaled))[seq_of_interest,])
  print((data %>% select(C26, C29, C31, {{C0_var}}, {{amp_var}}, {{phase_var}}))[seq_of_interest,])
  
  # phase_of_interest = 0
  # amp_of_interest = .5
  
  # phase_of_interest = 3*pi/8
  # amp_of_interest = .25
  
  sin_position_sequence = seq(min(position_sequence), max(position_sequence), by=0.5)
  # sin_position_sequence_calc = seq(from=-.5, length.out=length(sin_position_sequence), by=0.5)
  sin_position_sequence_calc = seq(from=-seq_of_interest*7+pi, length.out=length(sin_position_sequence), by=0.5)
  
  sine_31 = C0_of_interest + amp_of_interest*aa[3]*sin(2*pi*(sin_position_sequence_calc+31+9)/10.4 + phase_of_interest)
  sine_29 = C0_of_interest + amp_of_interest*aa[2]*sin(2*pi*(sin_position_sequence_calc+29+9)/10.4 + phase_of_interest)
  sine_26 = C0_of_interest + amp_of_interest*aa[1]*sin(2*pi*(sin_position_sequence_calc+26+9)/10.4 + phase_of_interest)
  
  alpha1=1
  
  ggplot(mapping=aes(x=position_sequence)) +
    # geom_line(aes(y=data$C26_rank, color="C26"), alpha=alpha1) +
    # geom_line(aes(y=data$C29_rank, color="C29"), alpha=alpha1) +
    # geom_line(aes(y=data$C31_rank, color="C31"), alpha=alpha1) +
    # geom_point(aes(y=data$C26_scaled, color="C26"), alpha=alpha1) +
    # geom_point(aes(y=data$C29_scaled, color="C29"), alpha=alpha1) +
    # geom_point(aes(y=data$C31_scaled, color="C31"), alpha=alpha1) +
    geom_point(aes(y=data$C26, color="C26"), alpha=alpha1) +
    geom_point(aes(y=data$C29, color="C29"), alpha=alpha1) +
    geom_point(aes(y=data$C31, color="C31"), alpha=alpha1) +
    geom_line(aes(x=sin_position_sequence, y=sine_31, color="C31"), alpha=0.3) + 
    geom_line(aes(x=sin_position_sequence, y=sine_29, color="C29"), alpha=0.3) + 
    geom_line(aes(x=sin_position_sequence, y=sine_26, color="C26"), alpha=0.3) + 
    xlab("Position") +
    ylab("Value") +
    scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green"))
}

random_center_point = sample(1:nrow(dat_tiling), size=1)
# plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=3)
# plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=5)
# plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=7)
# 
# plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=3,
#                                     C0_var=C0_lag1, phase_var=Phase_lag1, amp_var=Amplitude_lag1)
# plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=5,
#                                     C0_var=C0_lag1, phase_var=Phase_lag1, amp_var=Amplitude_lag1)
# plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=7,
#                                     C0_var=C0_lag1, phase_var=Phase_lag1, amp_var=Amplitude_lag1)

grid.arrange(plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=3),
             plot_rankings_C_by_library_position(dat_tiling, center_point = random_center_point, num_sequences=3, aa=c(1,1,1),
                                                 C0_var=C0_lag1, phase_var=Phase_lag1, amp_var=Amplitude_lag1))

library(gridExtra)
random_center_point = sample(1:nrow(chrv), size=1)
# plot_C_by_library_position(chrv, starting_point = random_starting_point, num_sequences=10, include_C0 = TRUE)
plot_rankings_C_by_library_position(chrv, center_point = random_center_point, num_sequences=3)
plot_rankings_C_by_library_position(chrv, center_point = random_center_point, num_sequences=5)
plot_rankings_C_by_library_position(chrv, center_point = random_center_point, num_sequences=7)


plot1_chrv_starting_point = 213500/7
plot1_chrv_noC0 = plot_C_by_library_position(chrv, starting_point = plot1_chrv_starting_point, num_sequences=50)
plot1_chrv_noC0
plot1_chrv_C0 = plot_C_by_library_position(chrv, starting_point = plot1_chrv_starting_point, num_sequences=50, include_C0=TRUE)
plot1_chrv_C0
grid.arrange(plot1_chrv_noC0, plot1_chrv_C0, nrow=2)

plot2_chrv_starting_point = 432100/7
plot2_chrv_noC0 = plot_C_by_library_position(chrv, starting_point = plot2_chrv_starting_point, num_sequences=50, include_C0 = FALSE)
plot2_chrv_C0 = plot_C_by_library_position(chrv, starting_point = plot2_chrv_starting_point, num_sequences=50, include_C0 = TRUE)
grid.arrange(plot2_chrv_noC0, plot2_chrv_C0, nrow=2)

# plot_C_by_library_position(chrv, starting_point = 547100/7, num_sequences=50)

# plot_C_by_library_position(tiling, starting_point = sample(1:nrow(tiling), size=1), num_sequences=50)
plot1_tiling_noC0 = plot_C_by_library_position(tiling, starting_point = 1, num_sequences=143)
plot1_tiling_noC0
plot1_tiling_C0 = plot_C_by_library_position(tiling, starting_point = 1, num_sequences=143, include_C0=TRUE)
plot1_tiling_C0
grid.arrange(plot1_tiling_noC0, plot1_tiling_C0, nrow=2)

# plot_C_by_library_position(tiling, starting_point = 143+1, num_sequences=143)
# plot_C_by_library_position(tiling, starting_point = 143*2+1, num_sequences=143)
# plot_C_by_library_position(tiling, starting_point = 143*575+1, num_sequences=143)

plot_C_by_library_position(tiling, starting_point = 143*sample(0:575, size=1)+1, num_sequences=143, include_C0=TRUE)



plot_smooth_function_of_individual_Cs = function(data, starting_point=1, num_sequences=10) {
  row_sequence = seq(from=starting_point, by=1, length.out=num_sequences)
  position_sequence = (row_sequence-1)*7
  # sin_position_sequence = seq((starting_point-1)*7, (starting_point-1)*7+7*(num_sequences-1))
  
  # temp_data = data.frame(index = position_sequence,
  #                        raw_Cvalues = c(data[row_sequence, "C26"], 
  #                                        data[row_sequence, "C29"],
  #                                        data[row_sequence, "C31"]),
  #                        group = rep(1:3, each = num_sequences))
  
  # df=num_sequences/3
  window_size=3
  
  C26_ma <- rollmean(data[row_sequence, "C26"], k = window_size, align = "center", 
                     fill = c(head(data[row_sequence, "C26"], 1), 0, tail(data[row_sequence, "C26"], 1)))
  C29_ma <- rollmean(data[row_sequence, "C29"], k = window_size, align = "center", 
                     fill = c(head(data[row_sequence, "C29"], 1), 0, tail(data[row_sequence, "C29"], 1)))
  C31_ma <- rollmean(data[row_sequence, "C31"], k = window_size, align = "center", 
                     fill = c(head(data[row_sequence, "C31"], 1), 0, tail(data[row_sequence, "C31"], 1)))
  
  # Plot with ggplot2 and geom_smooth
  # ggplot(temp_data, aes(x = index, y = raw_Cvalues, color = factor(group))) +
  #   geom_line(alpha=0.3) +
  #   # geom_smooth(method = "gam", formula = y ~ s(x, k=num_sequences)) +  # Use GAM smoothing
  #   # geom_smooth(method = "loess") +
  #   geom_smooth(method = "lm", formula = y ~ splines::bs(x, df = df),
  #               se=FALSE) +  
  #   labs(title = "Smoothed Function with Multiple Measurements",
  #        x = "Index", y = "Measurement")
  
  ggplot(mapping=aes(x=position_sequence, y=data[row_sequence,])) +
    geom_line(aes(y=data[row_sequence, "C26"], color="C26"), alpha=0.3) +
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C26"], color="C26"),
    #             method = "gam", formula = y ~ s(x, k=num_sequences)) +  # Use GAM smoothing
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C26"], color="C26"),
    #             method = "loess") +
    geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C26"], color="C26"),
                method = "lm", formula = y ~ splines::bs(x, df=df), se=FALSE) +
    geom_line(aes(y=data[row_sequence, "C29"], color="C29"), alpha=0.3) +
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C29"], color="C29"),
    #             method = "gam", formula = y ~ s(x, k=num_sequences)) +  # Use GAM smoothing
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C29"], color="C29"),
    #             method = "loess") +
    geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C29"], color="C29"),
                method = "lm", formula = y ~ splines::bs(x, df=df), se=FALSE) +
    geom_line(aes(y=data[row_sequence, "C31"], color="C31"), alpha=0.3) +
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C31"], color="C31"),
    #             method = "gam", formula = y ~ s(x, k=num_sequences)) +  # Use GAM smoothing
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C31"], color="C31"),
    #             method = "loess") +
    geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C31"], color="C31"),
                method = "lm", formula = y ~ splines::bs(x, df=df), se=FALSE) +
    labs(title = "Smoothed Individual C Values",
         x = "Position", y = "Value") +
    scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green"))
}

library(tidyquant)
plot_smooth_function_of_C0 = function(data, starting_point=1, num_sequences=10) {
  row_sequence = seq(from=starting_point, by=1, length.out=num_sequences)
  position_sequence = (row_sequence-1)*7
  
  # df=num_sequences/3
  window_size=3
  
  C0_ma <- rollmean(data[row_sequence, "C0"], k = window_size, align = "center", 
                    fill = c(head(data[row_sequence, "C0"], 1), 0, tail(data[row_sequence, "C0"], 1)))
  
  ggplot(mapping=aes(x=position_sequence, y=data[row_sequence,])) +
    geom_line(aes(y=data[row_sequence, "C26"], color="C26"), alpha=0.3) +
    geom_line(aes(y=data[row_sequence, "C29"], color="C29"), alpha=0.3) +
    geom_line(aes(y=data[row_sequence, "C31"], color="C31"), alpha=0.3) +
    # geom_smooth(mapping=aes(x=position_sequence, y=data[row_sequence, "C0"], color="C0"),
    #             method = "lm", formula = y ~ splines::bs(x, df=df), se=FALSE) +
    geom_line(aes(x=position_sequence, y=C0_ma, color="C0")) +
    labs(title = "Smoothed C0",
         x = "Position", y = "Value") +
    scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green", "C0"="brown"))
}

plot1_chrv_smooth = plot_smooth_function_of_individual_Cs(chrv, starting_point = plot1_chrv_starting_point, num_sequences=50)
grid.arrange(plot1_chrv_noC0, plot1_chrv_smooth, nrow=2)

plot1_chrv_smooth_C0 = plot_smooth_function_of_C0(chrv, starting_point = plot1_chrv_starting_point, num_sequences=50)
grid.arrange(plot1_chrv_noC0, plot1_chrv_smooth_C0, nrow=2)

grid.arrange(plot1_chrv_C0, plot1_chrv_smooth_C0, nrow=2)

plot2_chrv_smooth = plot_smooth_function_of_individual_Cs(chrv, starting_point = plot2_chrv_starting_point, num_sequences=50)
grid.arrange(plot2_chrv_noC0, plot2_chrv_smooth, nrow=2)
grid.arrange(plot2_chrv_noC0, plot2_chrv_C0, plot2_chrv_smooth, nrow=3)




plot1_tiling_smooth = plot_smooth_function_of_individual_Cs(tiling, starting_point = 1, num_sequences=143)
grid.arrange(plot1_tiling_noC0, plot1_tiling_smooth, nrow=2)

plot1_tiling_smooth_C0 = plot_smooth_function_of_C0(tiling, starting_point = 1, num_sequences=143)
grid.arrange(plot1_tiling_noC0, plot1_tiling_smooth_C0, nrow=2)

grid.arrange(plot1_tiling_C0, plot1_tiling_smooth_C0, nrow=2)

tiling_starting_point = 143*sample(0:575, size=1)+1
grid.arrange(plot_C_by_library_position(tiling, starting_point = tiling_starting_point, num_sequences=143, include_C0=TRUE),
             plot_smooth_function_of_C0(tiling, starting_point = tiling_starting_point, num_sequences=143))


grid.arrange(plot_smooth_function_of_C0(tiling, starting_point = tiling_starting_point, num_sequences=143),
             plot_smooth_function_of_C0(tiling, starting_point = tiling_starting_point, num_sequences=100))
