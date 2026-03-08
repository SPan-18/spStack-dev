# Calculate distance matrix

Computes the inter-site Euclidean distance matrix for one or two sets of
points.

## Usage

``` r
iDist(coords.1, coords.2, ...)
```

## Arguments

- coords.1:

  an \\n\times p\\ matrix with each row corresponding to a point in
  \\p\\-dimensional space.

- coords.2:

  an \\m\times p\\ matrix with each row corresponding to a point in
  \\p\\ dimensional space. If this is missing then `coords.1` is used.

- ...:

  currently no additional arguments.

## Value

The \\n\times n\\ or \\n\times m\\ inter-site Euclidean distance matrix.

## Author

Soumyakanti Pan <span18@ucla.edu>,  
Sudipto Banerjee <sudipto@ucla.edu>

## Examples

``` r
n <- 10
p1 <- cbind(runif(n),runif(n))
m <- 5
p2 <- cbind(runif(m),runif(m))
D <- iDist(p1, p2)
```
