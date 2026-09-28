# Round a stochastic vector or a row-stochastic matrix

Rounds a vector such that the sum of 1 is preserved. Rounds a matrix
such that each row sum up to 1. One entry is adjusted after rounding
such that the rounding error is the smallest.

## Usage

``` r
round_stochastic(x, digits = 7)
```

## Arguments

- x:

  a stochastic vector or a row-stochastic matrix.

- digits:

  number of digits for rounding.

## Value

The rounded vector or matrix.

## See also

[round](https://rdrr.io/r/base/Round.html)

Other utilities:
[`colors`](http://michael.hahsler.net/markovDP/reference/colors.md)

## Examples

``` r
# regular rounding would not sum up to 1
x <- c(0.333, 0.334, 0.333)

round_stochastic(x)
#> [1] 0.333 0.334 0.333
round_stochastic(x, digits = 2)
#> [1] 0.33 0.33 0.34
round_stochastic(x, digits = 1)
#> [1] 0.3 0.3 0.4
round_stochastic(x, digits = 0)
#> [1] 0 0 1


# round a stochastic matrix
m <- matrix(runif(15), ncol = 3)
m <- sweep(m, 1, rowSums(m), "/")

m
#>            [,1]      [,2]      [,3]
#> [1,] 0.41725134 0.3116414 0.2711072
#> [2,] 0.29966033 0.3461241 0.3542156
#> [3,] 0.29152543 0.4954197 0.2130549
#> [4,] 0.13446475 0.4366508 0.4288844
#> [5,] 0.09807003 0.6475371 0.2543929
round_stochastic(m, digits = 2)
#>      [,1] [,2] [,3]
#> [1,] 0.42 0.31 0.27
#> [2,] 0.30 0.35 0.35
#> [3,] 0.29 0.50 0.21
#> [4,] 0.13 0.44 0.43
#> [5,] 0.10 0.65 0.25
round_stochastic(m, digits = 1)
#>      [,1] [,2] [,3]
#> [1,]  0.4  0.3  0.3
#> [2,]  0.3  0.3  0.4
#> [3,]  0.3  0.5  0.2
#> [4,]  0.1  0.5  0.4
#> [5,]  0.1  0.6  0.3
round_stochastic(m, digits = 0)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    1    0    0
#> [3,]    1    0    0
#> [4,]    0    0    1
#> [5,]    0    1    0
```
