# Make two side-by-side surface plots

Make two side-by-side surface plots, particularly useful towards a
comparative study of two spatial surfaces.

## Usage

``` r
surfaceplot2(
  tab,
  coords_name,
  var1_name,
  var2_name,
  h = 8,
  col.pal,
  mark_points = FALSE
)
```

## Arguments

- tab:

  a data-frame containing spatial co-ordinates and the variables to plot

- coords_name:

  name of the two columns that contains the co-ordinates of the points

- var1_name:

  name of the column containing the first variable to be plotted

- var2_name:

  name of the column containing the second variable to be plotted

- h:

  integer; (optional) controls smoothness of the spatial interpolation
  as appearing in the
  [`MBA::mba.surf()`](https://rdrr.io/pkg/MBA/man/mba.surf.html)
  function. Default is 8.

- col.pal:

  Optional; color palette, preferably divergent, use `colorRampPalette`
  function from `grDevices`. Default is 'RdYlBu'.

- mark_points:

  Logical; if `TRUE`, the input points are marked. Default is `FALSE`.

## Value

a list containing two `ggplot` objects

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
data(simGaussian)
plots_2 <- surfaceplot2(simGaussian, coords_name = c("s1", "s2"),
                        var1_name = "z_true", var2_name = "y")
plots_2
#> [[1]]

#> 
#> [[2]]

#> 
```
