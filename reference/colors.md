# Default Colors for Visualization

Default discrete and continuous colors used in the package `markovDP`.

## Usage

``` r
colors_discrete(n, col = NULL)

colors_continuous(val, col = NULL)
```

## Arguments

- n:

  number of states.

- col:

  custom color palette. `colors_discrete()` uses the first n colors.
  `colors_continuous()` uses the given colors to calculate a palette
  (see
  [`grDevices::colorRamp()`](https://rdrr.io/r/grDevices/colorRamp.html)).
  The default is a blue-red color ramp.

- val:

  a vector with values to be translated to colors.

## Value

`colors_discrete()` returns a color palette and `colors_continuous()`
returns the colors associated with the supplied values.

## See also

Other utilities:
[`round_stochastic()`](http://michael.hahsler.net/markovDP/reference/round_stochastic.md)

## Examples

``` r
colors_discrete(5)
#> [1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00"

colors_continuous(runif(10))
#>  [1] "#E41A1C" "#8B6F8E" "#A56475" "#BA575E" "#D83632" "#377EB8" "#BA575D"
#>  [8] "#E31C1D" "#CD4644" "#837194"
```
