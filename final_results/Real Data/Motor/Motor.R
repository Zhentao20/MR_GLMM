setwd("~/Desktop/final_results/Motor")
load("simu220.08.RData")

library(fields)
library(RColorBrewer)
b_stimu <- bstimu[[1]]

mild_colors <- brewer.pal(9, "Blues")
par(mar = c(5, 4, 4, 8))  # Increase the right margin for the legend

image(t(b_stimu[nrow(b_stimu):1, ]), axes = FALSE, col = mild_colors)
axis(1, at = seq(0, 1, length.out = 50), labels = 1:50, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 50), labels = 50:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(b_stimu[nrow(b_stimu):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

b_stimu_f <- (b_stimu + t(b_stimu)) / 2
write.csv(b_stimu_f, "hcp_tongue.csv", row.names = T)
load("simu320.08.RData")

library(fields)
library(RColorBrewer)
the_stimu <- estheta

mild_colors <- brewer.pal(9, "Blues")
par(mar = c(5, 4, 4, 8))  # Increase the right margin for the legend

image(t(the_stimu[nrow(the_stimu):1, ]), axes = FALSE,  col = mild_colors)
axis(1, at = seq(0, 1, length.out = 50), labels = 1:50, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 50), labels = 50:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(the_stimu[nrow(the_stimu):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

write.csv(the_stimu, "hcp_intercept.csv", row.names = T)



