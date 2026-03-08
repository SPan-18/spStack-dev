# Create a collection of candidate models for stacking

Creates an object of class `'candidateModels'` that contains a list of
candidate models for stacking. The function takes a list of candidate
values for each model parameter and returns a list of possible
combinations of these values based on either simple aggregation or
Cartesian product of indivdual candidate values.

## Usage

``` r
candidateModels(params_list, aggregation = "simple")
```

## Arguments

- params_list:

  a list of candidate values for each model parameter. See examples for
  details.

- aggregation:

  a character string specifying the type of aggregation to be used.
  Options are `'simple'` and `'cartesian'`. Default is `'simple'`.

## Value

an object of class `'candidateModels'`

## See also

[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
m1 <- candidateModels(list(phi_s = c(1, 1), phi_t = c(1, 2)), "simple")
m1
#> [[1]]
#> [[1]]$phi_s
#> [1] 1
#> 
#> [[1]]$phi_t
#> [1] 1
#> 
#> 
#> [[2]]
#> [[2]]$phi_s
#> [1] 1
#> 
#> [[2]]$phi_t
#> [1] 2
#> 
#> 
#> attr(,"class")
#> [1] "candidateModels"
m2 <- candidateModels(list(phi_s = c(1, 1), phi_t = c(1, 2)), "cartesian")
m2
#> [[1]]
#> [[1]]$phi_s
#> [1] 1
#> 
#> [[1]]$phi_t
#> [1] 1
#> 
#> 
#> [[2]]
#> [[2]]$phi_s
#> [1] 1
#> 
#> [[2]]$phi_t
#> [1] 1
#> 
#> 
#> [[3]]
#> [[3]]$phi_s
#> [1] 1
#> 
#> [[3]]$phi_t
#> [1] 2
#> 
#> 
#> [[4]]
#> [[4]]$phi_s
#> [1] 1
#> 
#> [[4]]$phi_t
#> [1] 2
#> 
#> 
#> attr(,"class")
#> [1] "candidateModels"
m3 <- candidateModels(list(phi_s = list(c(1, 1), c(1, 2)),
                          phi_t = list(c(1, 3), c(2, 3)),
                          boundary = c(0.5, 0.75)),
                      "simple")
```
