setwd("~/Desktop/final_results/AAL")
load("simu220.1.RData")

library(fields)
library(RColorBrewer)
b1 <- bstimu[[1]]
b2 <- bstimu[[2]]
b3 <- bstimu[[3]]
b4 <- bstimu[[4]]

mild_colors <- brewer.pal(9, "Blues")
par(mar = c(5, 4, 4, 8))  # Increase the right margin for the legend

#b1
image(t(b1[nrow(b1):1, ]), axes = FALSE, col = mild_colors)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(b1[nrow(b1):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

#b2
image(t(b2[nrow(b2):1, ]), axes = FALSE, col = mild_colors)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(b2[nrow(b2):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

#b3
image(t(b3[nrow(b3):1, ]), axes = FALSE, col = mild_colors)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(b3[nrow(b3):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

#b4
image(t(b4[nrow(b4):1, ]), axes = FALSE, col = mild_colors)
axis(1, at = seq(0, 1, length.out = 50), labels = 1:50, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 50), labels = 50:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(b4[nrow(b4):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))


# Set up a 2x2 layout for four plots
par(mfrow = c(2, 2), mar = c(5, 4, 4, 6))  # Adjust margins for better spacing

# Plot 1: b1
image(t(b1[nrow(b1):1, ]), axes = FALSE, col = mild_colors, main = "Plot 1: Sex")
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)
image.plot(t(b1[nrow(b1):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

# Plot 2: b2
image(t(b2[nrow(b2):1, ]), axes = FALSE, col = mild_colors, main = "Plot 2: Age")
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)
image.plot(t(b2[nrow(b2):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

# Plot 3: b3
image(t(b3[nrow(b3):1, ]), axes = FALSE, col = mild_colors, main = "Plot 3: APOE4")
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)
image.plot(t(b3[nrow(b3):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

# Plot 4: b4
image(t(b4[nrow(b4):1, ]), axes = FALSE, col = mild_colors, main = "Plot 4: DX Label")
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)
image.plot(t(b4[nrow(b4):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))

# Reset layout
par(mfrow = c(1, 1))  # Reset to default layout

max(b1)
#1101
library(fields)  # Required for image.plot()

# Determine the global color scale limits based on all data
all_data <- c(b1, b2, b3, b4)
zlim <- range(all_data, na.rm = TRUE)

# Set up layout and margins
par(mfrow = c(2, 2), mar = c(5, 4, 4, 6))  # Adjust margins for better spacing

# Plot 1: b1
image(t(b1[nrow(b1):1, ]), axes = FALSE, col = mild_colors, main = "Plot 1: Sex", zlim = zlim)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Plot 2: b2
image(t(b2[nrow(b2):1, ]), axes = FALSE, col = mild_colors, main = "Plot 2: Age", zlim = zlim)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Plot 3: b3
image(t(b3[nrow(b3):1, ]), axes = FALSE, col = mild_colors, main = "Plot 3: APOE4", zlim = zlim)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Plot 4: b4
image(t(b4[nrow(b4):1, ]), axes = FALSE, col = mild_colors, main = "Plot 4: DX Label", zlim = zlim)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Shared legend for all plots
image.plot(legend.only = TRUE, col = mild_colors, zlim = zlim, legend.shrink = 0.8, 
           legend.width = 1.0, axis.args = list(cex.axis = 0.8))

# Reset layout
par(mfrow = c(1, 1))  # Reset to default layout

load("simu320.1.RData")

library(fields)
library(RColorBrewer)
the_stimu <- estheta

mild_colors <- brewer.pal(9, "Blues")
par(mar = c(5, 4, 4, 8))  # Increase the right margin for the legend

image(t(the_stimu[nrow(the_stimu):1, ]), axes = FALSE,  col = mild_colors)
axis(1, at = seq(0, 1, length.out = 90), labels = 1:90, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = 90), labels = 90:1, las = 2, cex.axis = 0.7)

# Add color legend using image.plot() from fields
image.plot(t(the_stimu[nrow(the_stimu):1, ]), add = TRUE, legend.only = TRUE, col = mild_colors,
           legend.shrink = 0.8, legend.width = 1.0, axis.args = list(cex.axis = 0.8))




