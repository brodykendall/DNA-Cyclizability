# Periodicity plots:
names <- c('Top Quartile', 'Bottom Quartile')
clrs <- c('red', 'blue')
ltype <- c(3, 4, 5, 1)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend = names, bty='n', fill=clrs, ncol=4, x.intersp=0.2)


# Unknown amplitude plots:
names <- c('C26', 'C29', 'C31', 'C0')
clrs <- c('blue', 'orange', 'green', 'red')
ltype <- c(3, 4, 5, 1)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend = names, bty='n', fill=clrs, ncol=4, x.intersp=0.2)


