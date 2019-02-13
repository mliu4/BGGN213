features <- read.delim("bimm143_05_rstats/feature_counts.txt")
par(mar=c(6.1,12.1,1.1,2.1))
barplot(features[,2], horiz = TRUE, xlim = c(0, 80000), xlab = "Feature count", names.arg = features[,1], main="Some Title", las = 1)