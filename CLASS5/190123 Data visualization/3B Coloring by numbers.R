genes <- read.table("bimm143_05_rstats/up_down_expression.txt", header = TRUE)
plot (genes[,2], genes[,3], col = genes$State, xlab = "Expression Condition 2", ylab = "Expression condition 3")
palette(c("red", "green", "blue"))
plot (genes[,2], genes[,3], col = genes$State, xlab = "Expression Condition 2", ylab = "Expression condition 3")