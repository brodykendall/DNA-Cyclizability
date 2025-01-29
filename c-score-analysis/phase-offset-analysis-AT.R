#### Top 4000/25000:

# Phase of fft at period=10 for A/T C26, in [-pi, pi]
Arg(find_ps1_fourier(C26, random %>% top_n(4000, C26)))
# 2.968705
Arg(find_ps1_fourier(C26, tiling %>% top_n(25000, C26)))
# 3.069499
Arg(find_ps1_fourier(C26, chrv %>% top_n(25000, C26)))
# -3.098992

# Phase of fft at period=10 for A/T C29, in [-pi, pi]
Arg(find_ps1_fourier(C29, random %>% top_n(4000, C29)))
# -1.244077
Arg(find_ps1_fourier(C29, tiling %>% top_n(25000, C29)))
# -1.373877
Arg(find_ps1_fourier(C29, chrv %>% top_n(25000, C29)))
# -1.364497

# Phase of fft at period=10 for A/T C31, in [-pi, pi]
Arg(find_ps1_fourier(C31, random %>% top_n(4000, C31)))
# -0.2807173
Arg(find_ps1_fourier(C31, tiling %>% top_n(25000, C31)))
# -0.4296983
Arg(find_ps1_fourier(C31, chrv %>% top_n(25000, C31)))
# -0.3889079


# Phase of fft at period=10 for A/T C26, in bp [0, 10]
(Arg(find_ps1_fourier(C26, random %>% top_n(4000, C26)))/(2*pi) + .5)*10
# 9.724841
(Arg(find_ps1_fourier(C26, tiling %>% top_n(25000, C26)))/(2*pi) + .5)*10
# 9.885259
(Arg(find_ps1_fourier(C26, chrv %>% top_n(25000, C26)))/(2*pi) + .5)*10
# 0.06780148

# Phase of fft at period=10 for A/T C29, in bp [0, 10]
(Arg(find_ps1_fourier(C29, random %>% top_n(4000, C29)))/(2*pi) + .5)*10
# 3.019989
(Arg(find_ps1_fourier(C29, tiling %>% top_n(25000, C29)))/(2*pi) + .5)*10
# 2.813406
(Arg(find_ps1_fourier(C29, chrv %>% top_n(25000, C29)))/(2*pi) + .5)*10
# 2.828336

# Phase of fft at period=10 for A/T C31, in bp [0, 10]
(Arg(find_ps1_fourier(C31, random %>% top_n(4000, C31)))/(2*pi) + .5)*10
# 4.553225
(Arg(find_ps1_fourier(C31, tiling %>% top_n(25000, C31)))/(2*pi) + .5)*10
# 4.316114
(Arg(find_ps1_fourier(C31, chrv %>% top_n(25000, C31)))/(2*pi) + .5)*10
# 4.381034





#### Top 2000/15000:

# Phase of fft at period=10 for A/T C26, in [-pi, pi]
Arg(find_ps1_fourier(C26, random %>% top_n(2000, C26)))
# 2.960222
Arg(find_ps1_fourier(C26, tiling %>% top_n(15000, C26)))
# 3.057243
Arg(find_ps1_fourier(C26, chrv %>% top_n(15000, C26)))
# -3.094566

# Phase of fft at period=10 for A/T C29, in [-pi, pi]
Arg(find_ps1_fourier(C29, random %>% top_n(2000, C29)))
# -1.26529
Arg(find_ps1_fourier(C29, tiling %>% top_n(15000, C29)))
# -1.369882
Arg(find_ps1_fourier(C29, chrv %>% top_n(15000, C29)))
# -1.346893

# Phase of fft at period=10 for A/T C31, in [-pi, pi]
Arg(find_ps1_fourier(C31, random %>% top_n(2000, C31)))
# -0.3296806
Arg(find_ps1_fourier(C31, tiling %>% top_n(15000, C31)))
# -0.4570286
Arg(find_ps1_fourier(C31, chrv %>% top_n(15000, C31)))
# -0.4126314


# Phase of fft at period=10 for A/T C26, in bp [0, 10]
(Arg(find_ps1_fourier(C26, random %>% top_n(2000, C26)))/(2*pi) + .5)*10
# 9.71134
(Arg(find_ps1_fourier(C26, tiling %>% top_n(15000, C26)))/(2*pi) + .5)*10
# 9.865754
(Arg(find_ps1_fourier(C26, chrv %>% top_n(15000, C26)))/(2*pi) + .5)*10
# 0.07484533

# Phase of fft at period=10 for A/T C29, in bp [0, 10]
(Arg(find_ps1_fourier(C29, random %>% top_n(2000, C29)))/(2*pi) + .5)*10
# 2.986229
(Arg(find_ps1_fourier(C29, tiling %>% top_n(15000, C29)))/(2*pi) + .5)*10
# 2.819765
(Arg(find_ps1_fourier(C29, chrv %>% top_n(15000, C29)))/(2*pi) + .5)*10
# 2.856353


# Phase of fft at period=10 for A/T C31, in bp [0, 10]
(Arg(find_ps1_fourier(C31, random %>% top_n(2000, C31)))/(2*pi) + .5)*10
# 4.475297
(Arg(find_ps1_fourier(C31, tiling %>% top_n(15000, C31)))/(2*pi) + .5)*10
# 4.272616
(Arg(find_ps1_fourier(C31, chrv %>% top_n(15000, C31)))/(2*pi) + .5)*10
# 4.343277


