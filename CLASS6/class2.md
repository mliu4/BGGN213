Class 6 R Functions
================
MENGDAN LIU
January 25, 2019

Time to learn about **read.table** and its friends.

Inserting a code chunk

``` r
read.table("https://bioboot.github.io/bggn213_S18/class-material/test1.txt", header = TRUE, sep = ",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file1 <- "https://bioboot.github.io/bggn213_S18/class-material/test1.txt"
read.csv(file1)
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file2 <- "https://bioboot.github.io/bggn213_S18/class-material/test2.txt"
file3 <- "https://bioboot.github.io/bggn213_S18/class-material/test3.txt"
read.table(file2, header = TRUE, sep = "$")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table(file3)
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

R functions
-----------

``` r
add <- function(x, y=1){
  return(x+y)
}
```

Let's use our useless function

``` r
#Return 2
add(1)
```

    ## [1] 2

``` r
#Return 3
add(1+2)
```

    ## [1] 4

``` r
#It can do vectors too
add(c(1,2,3,4))
```

    ## [1] 2 3 4 5

``` r
#WOW
add(c(1,2,3,4), c(1,2,3,4))
```

    ## [1] 2 4 6 8

Why use a function?
===================

Functions let us repeat important operations without having to type them out every time. It's like getting a cloth towel to wipe up spills instead of using a whole bunch of paper towel rolls. YOu can always wash the cloth towel and use it again for the same purpose.

Here's an example function

``` r
Normalize <- function(inputDf, col){
  inputDf <- (inputDf$col-min(inputDf$col))/(max(inputDf$col)-min(inputDf$col))
  return(inputDf)
}
```

But this can be done better. After some optimizations, we get this:

``` r
rescale2 <- function(x){
  rng <- range(x)
  x <- (x- rng[1])/(rng[2]-rng[1])
  return <- x
}
```

Further optimizations can be acconted for. Like covering NA data entries

``` r
rescale2 <- function(x){
  rng <- range(x, na.rm = TRUE)
  x <- (x-rng[1])/(rng[2]-rng[1])
  return(x)
}
```

We can do even better.

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
 if(na.rm) {
 rng <-range(x, na.rm=na.rm)
 } else {
 rng <-range(x)
 }
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 return(answer)
}
```

Okay, time for in class lab stuff.

HANDS ON LAB SECTION
--------------------

SECTION 1A
----------

``` r
df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)

fixer <- function(x, na.rm=TRUE){
  rnger <- range(x)
  x <- (x-rnger[1])/(rnger[2]-rnger[1])
  return(x)
}

for(i in df){
  i <- fixer(i)
  print(i)
}
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000
    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000
    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000
    ##  [1] NA NA NA NA NA NA NA NA NA NA

SECTION 1B
----------

``` r
#Function relies on functions from bio3d library.
library(bio3d)

#plotBFactor is a function that accepts a PDB iD as argument and outputs a graph of the corresponding protein's B-factor or the amount of wiggle atoms have in a structure(?)
plot.Achain.BFactor <- function(pdbNum){
  #Assigning pdb coordinate file to variable s
  inPDB <- read.pdb(pdbNum)
  #Assigns chainA atom coordinates to var s.chainA using trim.pdb
  inPDB.chainA <- trim.pdb(inPDB, chain="A", elety="CA")
  #Assigning chainA atom b values to s.b
  inPDB.bVals <- inPDB.chainA$atom$b
  #Plotting b values against chainA secondary structure using plot3b
  plotb3(inPDB.bVals, sse = inPDB.chainA, typ="l", ylab="Bfactor")
  
}
plot.Achain.BFactor("4AKE")
```

    ##   Note: Accessing on-line PDB file

![](class2_files/figure-markdown_github/unnamed-chunk-11-1.png)
