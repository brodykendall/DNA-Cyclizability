library(tidyverse)

library(ggplot2)
library(patchwork)

# dat_chrI = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrI_1bpresolution_subsequence50_predictions.csv")
# dat_chrII = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrII_1bpresolution_subsequence50_predictions.csv")
# dat_chrIII = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrIII_1bpresolution_subsequence50_predictions.csv")
# dat_chrIV = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrIV_1bpresolution_subsequence50_predictions.csv")
dat_chrV = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_predictions.csv")
# dat_chrVI = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrVI_1bpresolution_subsequence50_predictions.csv")
# dat_tiling = read.csv("data/Created/ir_lstm_cn_tiling_tiling_full_reconstructed_predictions.csv")

# N = nrow(dat_chrI)
# N = nrow(dat_chrII)
# N = nrow(dat_chrIII)
# N = nrow(dat_chrIV)
N = nrow(dat_chrV)
# N = nrow(dat_chrVI)
## N = nrow(dat_tiling)
# N=1005

# f26=fft(dat_chrI$n.26)
# f29=fft(dat_chrI$n.29)
# f31=fft(dat_chrI$n.31)

# f26=fft(dat_chrII$n.26)
# f29=fft(dat_chrII$n.29)
# f31=fft(dat_chrII$n.31)

# f26=fft(dat_chrIII$n.26)
# f29=fft(dat_chrIII$n.29)
# f31=fft(dat_chrIII$n.31)

# f26=fft(dat_chrIV$n.26)
# f29=fft(dat_chrIV$n.29)
# f31=fft(dat_chrIV$n.31)

f26=fft(dat_chrV$n.26)
f29=fft(dat_chrV$n.29)
f31=fft(dat_chrV$n.31)
# f26=fft(dat_chrV$n.26_rev_comp)
# f29=fft(dat_chrV$n.29_rev_comp)
# f31=fft(dat_chrV$n.31_rev_comp)

# f26=fft(dat_chrVI$n.26)
# f29=fft(dat_chrVI$n.29)
# f31=fft(dat_chrVI$n.31)

# tiling_gene_number = 7
# f26=fft(dat_tiling$n.26[((tiling_gene_number-1)*1005 + 1):(tiling_gene_number*1005)])
# f29=fft(dat_tiling$n.29[((tiling_gene_number-1)*1005 + 1):(tiling_gene_number*1005)])
# f31=fft(dat_tiling$n.31[((tiling_gene_number-1)*1005 + 1):(tiling_gene_number*1005)])

freqs = ((1:floor(N/2))-1)/(N)



fourier_df = data.frame("Magnitude_C26" = Mod(f26[1:floor(N/2)]),
                        "Magnitude_C29" = Mod(f29[1:floor(N/2)]),
                        "Magnitude_C31" = Mod(f31[1:floor(N/2)]),
                        "Frequency" = freqs)



# alpha_val = 0.1
alpha_val = 0.3
ylims=c(0,10000)
line_width=0.5
C26_freq_plot = ggplot(data=fourier_df, aes(x=Frequency)) +
  geom_line(aes(y=Magnitude_C26), color="blue", alpha=alpha_val, linewidth=line_width) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(n.breaks=3, limits=ylims)
C29_freq_plot = ggplot(data=fourier_df, aes(x=Frequency)) +
  geom_line(aes(y=Magnitude_C29), color="orange", alpha=alpha_val, linewidth=line_width) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(n.breaks=3, limits=ylims)
C31_freq_plot = ggplot(data=fourier_df, aes(x=Frequency)) +
  geom_line(aes(y=Magnitude_C31), color="green", alpha=alpha_val, linewidth=line_width) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(n.breaks=3, limits=ylims)



# alpha_val = 0.1
# C26_freq_plot = ggplot(data=fourier_df, aes(x=Frequency)) +
#   geom_line(aes(y=Amplitude_C26), color="blue", alpha=alpha_val)
# C29_freq_plot = ggplot(data=fourier_df, aes(x=Frequency)) +
#   geom_point(aes(y=Amplitude_C29), color="orange", alpha=alpha_val)
# C31_freq_plot = ggplot(data=fourier_df, aes(x=Frequency)) +
#   geom_point(aes(y=Amplitude_C31), color="green", alpha=alpha_val)

combined_freq_plot = C26_freq_plot / C29_freq_plot / C31_freq_plot
combined_freq_plot

# # find maximum amplitude of the fourier transform
# m26<- which.max(abs(f26))
# m29<- which.max(abs(f29))
# m31<- which.max(abs(f31))
# 
# #find phase corresponding to the maximum Fourier transform
# p26= atan2(Im(f26[m26]), Re(f26[m26]))
# p29= atan2(Im(f26[m29]), Re(f26[m29]))
# p31= atan2(Im(f26[m31]), Re(f26[m31]))
# 
# k26=2*pi *(29-26)/(p29-p26)
# print(k26)
# k29=2*pi *(31-26)/(p31-p26)
# print(k29)
# k31=2*pi *(31-29)/(p31-p29)
# print(k31)
# 
# N/(m26-1)
# N/(m29-1)
# N/(m31-1)



# find maximum amplitude of the fourier transform (for frequency > 0.05)
valid_indices = which(freqs > 0.05)

m26_valid = which.max(abs(f26)[valid_indices])
m29_valid = which.max(abs(f29)[valid_indices])
m31_valid = which.max(abs(f31)[valid_indices])

#find phase corresponding to the maximum (valid) Fourier transform
p26_valid = atan2(Im(f26[valid_indices][m26_valid]), Re(f26[valid_indices][m26_valid]))
p29_valid = atan2(Im(f29[valid_indices][m29_valid]), Re(f29[valid_indices][m29_valid]))
p31_valid = atan2(Im(f31[valid_indices][m31_valid]), Re(f31[valid_indices][m31_valid]))

# 1/freqs[valid_indices[m26_valid]]
# 1/freqs[valid_indices[m29_valid]]
# 1/freqs[valid_indices[m31_valid]]

k26_29_valid=2*pi *(26-29)/((p29_valid-p26_valid)%%(2*pi) - 2*pi)
print(k26_29_valid)
# I: 10.2669
# II: 10.57223
# III: 10.80559
# IV: 10.61351
# V: 11.06557
# VI: 5.143499
k26_31_valid=2*pi *(26-31)/((p31_valid-p26_valid)%%(2*pi) - 2*pi)
print(k26_31_valid)
# I: 10.98841
# II: 11.29074
# III: 6.85253
# IV: 11.49581
# V: 11.84209
# VI: 6.794949
k29_31_valid=2*pi *(29-31)/((p31_valid-p29_valid)%%(2*pi) - 2*pi)
print(k29_31_valid)
# I: 12.28323
# II: 12.5724
# III: 4.424549
# IV: 13.13351
# V: 13.23526
# VI: 13.10787



N/(valid_indices[m26_valid]-1)
# I: 10.31408
# II: 10.25326
# III: 10.23442
# IV: 10.22278
# V: 10.596
# VI: 10.30922
N/(valid_indices[m29_valid]-1)
# I: 10.31408
# II: 10.25326
# III: 10.23442
# IV: 10.22278
# V: 10.596
# VI: 10.11163
N/(valid_indices[m31_valid]-1)
# I: 10.31408
# II: 10.25326
# III: 10.32218
# IV: 10.22278
# V: 10.596
# VI: 10.11163

tiling_gene_fourier = lapply(1:576, function(gene_number) {
  f26=fft(dat_tiling$n.26[((gene_number-1)*1005 + 1):(gene_number*1005)])
  f29=fft(dat_tiling$n.29[((gene_number-1)*1005 + 1):(gene_number*1005)])
  f31=fft(dat_tiling$n.31[((gene_number-1)*1005 + 1):(gene_number*1005)])
  
  m26_valid = which.max(abs(f26)[valid_indices])
  m29_valid = which.max(abs(f29)[valid_indices])
  m31_valid = which.max(abs(f31)[valid_indices])
  
  ret = list()
  
  ret[[1]] = m26_valid
  ret[[2]] = m29_valid
  ret[[3]] = m31_valid
  
  ret[[4]] = N/(valid_indices[m26_valid]-1)
  ret[[5]] = N/(valid_indices[m29_valid]-1)
  ret[[6]] = N/(valid_indices[m31_valid]-1)
  
  ret[[7]] = f26
  ret[[8]] = f29
  ret[[9]] = f31
  
  return(ret)
})

freq_26_list = map(tiling_gene_fourier, ~freqs[valid_indices][.x[[1]]])
freq_29_list = map(tiling_gene_fourier, ~freqs[valid_indices][.x[[2]]])
freq_31_list = map(tiling_gene_fourier, ~freqs[valid_indices][.x[[3]]])

period_26_list = map(tiling_gene_fourier, ~.x[[4]])
period_29_list = map(tiling_gene_fourier, ~.x[[5]])
period_31_list = map(tiling_gene_fourier, ~.x[[6]])

mean(unlist(freq_26_list))
# 0.09527709
mean(unlist(freq_29_list))
# 0.09586616
mean(unlist(freq_31_list))
# 0.09532891

mean(c(mean(unlist(freq_26_list)), mean(unlist(freq_29_list)), mean(unlist(freq_31_list))))
# 0.09549072
1/mean(c(mean(unlist(freq_26_list)), mean(unlist(freq_29_list)), mean(unlist(freq_31_list))))
# 10.47222

1/mean(unlist(freq_26_list))
# 10.4957
1/mean(unlist(freq_29_list))
# 10.43121
1/median(unlist(freq_31_list))
# 10.46875

mean(c(1/mean(unlist(freq_26_list)), 1/mean(unlist(freq_29_list)), 1/median(unlist(freq_31_list))))
# 10.46522

mean(unlist(period_26_list))
# 10.52264
mean(unlist(period_29_list))
# 10.45608
mean(unlist(period_31_list))
# 10.51481


phase_diff_29_26_list = map(tiling_gene_fourier, ~c((Arg(.x[[8]][valid_indices][.x[[2]]]) - 
                                                       Arg(.x[[7]][valid_indices][.x[[1]]]))))
phase_diff_31_26_list = map(tiling_gene_fourier, ~c((Arg(.x[[9]][valid_indices][.x[[3]]]) - 
                                                       Arg(.x[[7]][valid_indices][.x[[1]]]))))
phase_diff_31_29_list = map(tiling_gene_fourier, ~c((Arg(.x[[9]][valid_indices][.x[[3]]]) - 
                                                       Arg(.x[[8]][valid_indices][.x[[2]]]))))


2*pi*(26-29)/median(unlist(phase_diff_29_26_list)%%(2*pi)-2*pi)
# 10.43183
2*pi*(26-31)/median(unlist(phase_diff_31_26_list)%%(2*pi)-2*pi)
# 11.30778
2*pi*(29-31)/median(unlist(phase_diff_31_29_list)%%(2*pi)-2*pi)
# 12.8736


C26_all_mags_df = do.call(cbind.data.frame, map(tiling_gene_fourier, ~Mod(.x[[7]][1:floor(N/2)])))
C29_all_mags_df = do.call(cbind.data.frame, map(tiling_gene_fourier, ~Mod(.x[[8]][1:floor(N/2)])))
C31_all_mags_df = do.call(cbind.data.frame, map(tiling_gene_fourier, ~Mod(.x[[9]][1:floor(N/2)])))





fourier_tiling_df = data.frame("Magnitude_C26" = rowMeans(C26_all_mags_df),
                               "Magnitude_C29" = rowMeans(C29_all_mags_df),
                               "Magnitude_C31" = rowMeans(C31_all_mags_df),
                               "Frequency" = freqs)

alpha_val = 1
ylims=c(0,130)
line_width=0.5
C26_freq_plot = ggplot(data=fourier_tiling_df, aes(x=Frequency)) +
  geom_line(aes(y=Magnitude_C26), color="blue", alpha=alpha_val, linewidth=line_width) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  ylim(ylims) +
  scale_y_continuous(n.breaks=2)
C29_freq_plot = ggplot(data=fourier_tiling_df, aes(x=Frequency)) +
  geom_line(aes(y=Magnitude_C29), color="orange", alpha=alpha_val, linewidth=line_width) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  ylim(ylims) +
  scale_y_continuous(n.breaks=2)
C31_freq_plot = ggplot(data=fourier_tiling_df, aes(x=Frequency)) +
  geom_line(aes(y=Magnitude_C31), color="green", alpha=alpha_val, linewidth=line_width) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ylim(ylims) +
  scale_y_continuous(n.breaks=2)

combined_freq_plot = C26_freq_plot / C29_freq_plot / C31_freq_plot
combined_freq_plot
