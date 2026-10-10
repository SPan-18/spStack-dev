# Make a surface plot

Make a surface plot

## Usage

``` r
surfaceplot(tab, coords_name, var_name, h = 8, col.pal, mark_points = FALSE)
```

## Arguments

- tab:

  a data-frame containing spatial co-ordinates and the variable to plot

- coords_name:

  name of the two columns that contains the co-ordinates of the points

- var_name:

  name of the column containing the variable to be plotted

- h:

  integer; (optional) controls smoothness of the spatial interpolation
  as appearing in the
  [`MBA::mba.surf()`](https://finleya.github.io/MBA/reference/mba.surf.html)
  function. Default is 8.

- col.pal:

  Optional; color palette, preferably divergent, use `colorRampPalette`
  function from `grDevices`. Default is the colorblind-friendly
  diverging palette 'RdBu' from ColorBrewer.

- mark_points:

  Logical; if `TRUE`, the input points are marked. Default is `FALSE`.

## Value

a `ggplot` object containing the surface plot

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
data(simSpatial)
plot1 <- surfaceplot(simSpatial, coords_name = c("s1", "s2"),
                     var_name = "z_true")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the spStack package.
#>   Please report the issue at <https://github.com/SPan-18/spStack-dev/issues>.
plot1


# try your favourite color palette
plot2 <- surfaceplot(simSpatial, coords_name = c("s1", "s2"),
                     var_name = "z_true",
                     col.pal = hcl.colors(100, "PuOr", rev = TRUE))
plot2
```
