Untitled
================

MXL genotype analysis
---------------------

How many G|G genotypes are there, in the MXL population

``` r
##Reading in mxl population genome file.
mxl <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

``` r
gg.pop <- mxl$Genotype..forward.strand.
table(gg.pop) / nrow(mxl) *100
```

    ## gg.pop
    ##     A|A     A|G     G|A     G|G 
    ## 34.3750 32.8125 18.7500 14.0625

FastQ quality scores
--------------------

``` r
#install.packages("seqinr")
#install.packages("gtools")

library("seqinr")
library("gtools")

chars <- s2c("DDDDCDEDCDDDDBBDDDCC@")
chars
```

    ##  [1] "D" "D" "D" "D" "C" "D" "E" "D" "C" "D" "D" "D" "D" "B" "B" "D" "D"
    ## [18] "D" "C" "C" "@"

``` r
phred <- asc(chars)-33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Working on the final file
-------------------------

``` r
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
table(expr$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
inds.gg <- expr$geno=="G/G"
summary(expr[inds.gg,"exp"])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   6.675  16.903  20.074  20.594  24.457  33.956

``` r
inds.ag <- expr$geno=="A/G"
summary(expr[inds.ag,"exp"])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
inds.aa <- expr$geno=="A/A"
summary(expr[inds.aa,"exp"])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   11.40   27.02   31.25   31.82   35.92   51.52

``` r
boxplot(exp~geno, data=expr)
```

![](CLASS13_files/figure-markdown_github/unnamed-chunk-6-1.png)
