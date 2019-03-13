CLASS7
================
Mengdan Liu
January 30, 2019

Functions revisited
-------------------

Load a source file with our rescale function from last time.

``` r
source("http://tinyurl.com/rescale-R")
```

Test this function

``` r
rescale(1:5)
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
#This will deliberately throw an error, for educational purposes and thus has been commented out
#rescale(c(1:5, "string"))
```

We want to make this function more robust to these types of errors.

``` r
print("I wish I had access to the updated version")
```

    ## [1] "I wish I had access to the updated version"

Let's try something else. Something with a new function.

``` r
#Let's define an example x and y
x <- c(1,2,NA,3,NA)
y <- c(NA,3,NA,3,4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
which(is.na(x)&is.na(y))
```

    ## [1] 3

Let's turn it into a function

``` r
bothNA <- function(x,y){
  sum(is.na(x)&is.na(y))
  
  
}
bothNA(x,y)
```

    ## [1] 1
