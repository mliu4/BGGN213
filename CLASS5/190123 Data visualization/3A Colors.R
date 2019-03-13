genderwar <-  "bimm143_05_rstats/male_female_counts.txt"
boysGirls <- read.delim(genderwar, header = TRUE)
barplot(boysGirls[,2], names.arg = boysGirls[,1], las = 2, col = rainbow(nrow(boysGirls)))
barplot(boysGirls[,2], names.arg=boysGirls[,1], las = 2, col= c("blue", "red"))
