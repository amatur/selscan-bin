data_rehh <- read.table("/home/art/workspace/selscan-bin/data/rehh/out_ihh_20k_rehh", header = FALSE)
data_rehh <- data.frame(
  x = data_rehh[, 1],             # X-axis values
  y1 = data_rehh[, 4],       # Y-axis values for the first line
  y2 = data_rehh[, 5]        # Y-axis values for the second line
)


data_hapbin <- read.table("/home/art/workspace/hapbin_ihh", header = FALSE)
data_hapbin <- data.frame(
  x = data_hapbin[, 1],             # X-axis values
  y1 = data_hapbin[, 4],       # Y-axis values for the first line
  y2 = data_hapbin[, 5]        # Y-axis values for the second line
)

#selbin++ nsl
data_nsl <- read.table("/home/art/workspace/selscan-bin/outbin.nsl.out", header = FALSE)
data_nsl <- data.frame(
  x = data_nsl[, 1],             # X-axis values
  y1 = data_nsl[, 4],       # Y-axis values for the first line
  y2 = data_nsl[, 5]        # Y-axis values for the second line
)

data_hapbin <- merge(data_nsl, data_hapbin, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE)
# print(xy)

# selscan-bin/data/out.20k.simple.map
data_selscan <- read.table("/home/art/workspace/selscan-bin/data/validation/out20k.ihs.out", header = FALSE)
#data_selscan <- read.table("/home/art/workspace/selscan-bin/data/outfile.nsl.out", header = FALSE)
data_selscan <- data.frame(
  x = data_selscan[, 1],             # X-axis values
  y1 = data_selscan[, 4],       # Y-axis values for the first line
  y2 = data_selscan[, 5]        # Y-axis values for the second line
)

#data2 <- read.table("/home/art/workspace/selscan-bin/data/validation/outbin20k.ihh", header = FALSE, sep = " ")


#data_ihh_curr <- read.table("/home/art/workspace/selscan-bin/selbin.ihh.distance1", header = FALSE, sep = " ")
data_ihh_curr <- read.table("/home/art/workspace/selscan-bin/selbin.ihh.distance_act", header = FALSE, sep = " ")
#data_ihh_curr <- read.table("/home/art/workspace/selscan-bin/outbin.ihh.out", header = FALSE, sep = " ")
data_ihh_curr <- data.frame(
  x = data_ihh_curr[, 1],             # X-axis values
  y1 = data_ihh_curr[, 4]/1,       # Y-axis values for the first line
  y2 = data_ihh_curr[, 5]/1        # Y-axis values for the second line
)

data2 <- read.table("/home/art/workspace/selscan-bin/outbin.ihh.out.long", header = FALSE, sep = " ")
data2 <- data.frame(
  x = data2[, 1],             # X-axis values
  y1 = data2[, 4],       # Y-axis values for the first line
  y2 = data2[, 5]        # Y-axis values for the second line
)

# Create a blank plot with appropriate axis limits
plot(NULL, xlim = c(5200, 5400), ylim = range(data_ihh_curr$y1, data_ihh_curr$y2*2), 
     xlab = "Locus", ylab = "iHH1", main = "Comparing iHH1 values (chr10: 20,000 sites, 6404 haps) ")
#retained 1502

lines(data_selscan$x, data_selscan$y1, col = "blue", lty = 1, lwd = 2)
lines(data_ihh_curr$x, data_ihh_curr$y1, col = "red", lty = 2, lwd = 2)
# lines(data_nsl$x, data_hapbin$y1.y, col = "magenta", lty = 3, lwd = 2)
# 

# xy <- merge(data_nsl, data_hapbin, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE)

# lines(data_nsl$x, data_nsl$y1, col = "black", lty = 3, lwd = 2)
lines(data_rehh$x, data_rehh$y1, col = "black", lty = 3, lwd = 2)


# Add a legend
legend("topleft", legend = c("Selscan iHH", "Selscan-NEW iHH", "Hapbin","Selscan-NEW nSL"), col = c("blue", "red", "magenta","black"), 
       lty = c(1, 2,3,4), cex=0.5,lwd = 2, title = "Tool name")
