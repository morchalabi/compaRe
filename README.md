# compaRe
To measure similarity between samples, we developed COMPARE. It measures the similarity between two datasets with any number of dimensions and observations.

## Installation

Run the following snippet in R or RStudio to install the package:

``` r
if(!require(devtools)
{
  install.packages("devtools")
}
if(!require(igraph)
{
  install.packages("igraph")
}
devtools::install_github(repo = 'morchalabi/compaRe', ref = 'dev', dependencies = T)
```
