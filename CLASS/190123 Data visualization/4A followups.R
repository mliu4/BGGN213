source("../bimm143_05_rstats/color_to_value_map.r")

meth <- read.delim("../bimm143_05_rstats/expression_methylation.txt", header = TRUE)
myColors <- map.colors(meth$expression,c(max(meth$expression), min(meth$expression)), colorRampPalette(c("blue", "red"))(100))

plot(meth$promoter.meth, meth$gene.meth, xlab = "Promoter Methylation", ylab = "Gene Methylation", col = myColors)

