# Create a sample dataset with two columns


data <- read.table("/home/art/workspace/pasted", header = FALSE, sep = " ")
data <- data.frame(
  x = 1:216,             # X-axis values
  y1 = data[, 1],       # Y-axis values for the first line
  y2 = data[, 2]        # Y-axis values for the second line
)

# Create a blank plot with appropriate axis limits
plot(NULL, xlim = c(1, 216), ylim = range(data$y1, data$y2), 
     xlab = "Decay of EHH from locus 0", ylab = "EHH (core haplotype 1)", main = "Checking that values match for chr10 20000x6404 ")

# Add the first line in blue
lines(data$x, data$y1, col = "blue", lty = 1, lwd = 2)

# Add the second line in red
lines(data$x, data$y2, col = "red", lty = 2, lwd = 2)

# Add a legend
legend("topright", legend = c("Selscan-NEW", "Selscan"), col = c("blue", "red"), 
       lty = c(1, 2), lwd = 2, title = "Tool name")
