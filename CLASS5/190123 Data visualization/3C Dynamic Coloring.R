meth = read.table("bimm143_05_rstats/expression_methylation.txt", header = TRUE)
color <- densCols(meth$gene.meth[indices], meth$expression[indices])
color.custom <- densCols(meth$gene.meth[indices], meth$expression[indices], colramp = colorRampPalette(c("blue2", "green2", "red2", "yellow")))
plot(meth$gene.meth[indices], meth$express[indices], col = color.custom, pch = 20)