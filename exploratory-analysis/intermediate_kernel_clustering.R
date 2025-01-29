library(cluster)
tiling_first_conv_kernels.conv_only.km = kmeans(tiling_first_conv_kernels.conv_only, 4)
clusplot(tiling_first_conv_kernels.conv_only,
         tiling_first_conv_kernels.conv_only.km$cluster,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         plotchar = FALSE,
         labels = 2,
         span = TRUE,
         main = paste("Clusters of first convolutional layer (3x4), 48 features"))



tiling_first_conv_kernels.conv_only_2.km = kmeans(tiling_first_conv_kernels.conv_only_2, 4)
clusplot(tiling_first_conv_kernels.conv_only_2,
         tiling_first_conv_kernels.conv_only_2.km$cluster,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         plotchar = FALSE,
         labels = 2,
         span = TRUE,
         main = paste("Clusters of first convolutional layer (3x4), 16 features"))
