library(cluster)
tiling_tiling_first_conv_output.conv_only.km = kmeans(tiling_tiling_first_conv_output.conv_only, 4)
clusplot(tiling_tiling_first_conv_output.conv_only,
         tiling_tiling_first_conv_output.conv_only.km$cluster,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         plotchar = FALSE,
         labels = 2,
         span = TRUE,
         main = paste("Clusters of output from first convolutional layer (3x4), 48 features"))
