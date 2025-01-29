ir_lstm_fits_random <- read.csv("~/Documents/Northwestern/Jiping/benchmarks/ir_lstm_fits_random.txt")
colnames(ir_lstm_fits_random) = c("accuracy", "MSE", "mean", "std", "reverse_accuracy")

# X1: yeast
# X3: random
# X5: tiling
# X6: chrV
ir_lstm_fits_random_yeast = ir_lstm_fits_random[(1:10)*4 - 3,]
ir_lstm_fits_random_val = ir_lstm_fits_random[(1:10)*4 - 2,]
ir_lstm_fits_random_tiling = ir_lstm_fits_random[(1:10)*4 - 1,]
ir_lstm_fits_random_chrV = ir_lstm_fits_random[(1:10)*4,]

max(ir_lstm_fits_random_yeast$accuracy)
max(ir_lstm_fits_random_val$accuracy)
max(ir_lstm_fits_random_tiling$accuracy)
max(ir_lstm_fits_random_chrV$accuracy)



ir_lstm_fits_tiling <- read.csv("~/Documents/Northwestern/Jiping/benchmarks/ir_lstm_fits_tiling.txt")
colnames(ir_lstm_fits_tiling) = c("accuracy", "MSE", "mean", "std", "reverse_accuracy")

ir_lstm_fits_tiling_yeast = ir_lstm_fits_tiling[(1:10)*4 - 3,]
ir_lstm_fits_tiling_random = ir_lstm_fits_tiling[(1:10)*4 - 2,]
ir_lstm_fits_tiling_val = ir_lstm_fits_tiling[(1:10)*4 - 1,]
ir_lstm_fits_tiling_chrV = ir_lstm_fits_tiling[(1:10)*4,]

max(as.numeric(ir_lstm_fits_tiling_yeast$accuracy))
max(as.numeric(ir_lstm_fits_tiling_random$accuracy))
max(as.numeric(ir_lstm_fits_tiling_val$accuracy))
max(as.numeric(ir_lstm_fits_tiling_chrV$accuracy))

