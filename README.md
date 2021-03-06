# COMPARE
COMPARE is an R package with two functions for similarity measurement and clustering. It employs a dynamic programming algorithm to measure spatial similarity between two given datasets. One can then use this function to generate a similarity/affinity matrix which can then be input to the clustering function. The clustering function reports maximal cliques as clusters and dense areas as communities. COMPARE
1. is data-size independent
1. can effectively circumvent background noise (like batch effect)
1. does not need  dimension reduction
1. does not need subsampling

Read the accompanying help page of the package for more details about setup and running.

For more information about its algorithms check out COMPARE-suite's [similarity](https://github.com/morchalabi/COMPARE-suite/wiki/COMPARE-Suite#similarity-matrix-generator) and [clustering](https://github.com/morchalabi/COMPARE-suite/wiki/COMPARE-Suite#clustering) modules. COMPARE-suite is a pipeline for high-throughput and high-content multiparameter screening built on COMPARE.

## Installation

Open the R interpreter or [RStudio](https://rstudio.com/products/rstudio/download/) IDE and type in the following command:

    if(!requireNamespace("igraph"))
    {
      message('<< Installing igraph >>')
      install.packages('igraph', repos = "https://cloud.r-project.org")
    }
    if(!requireNamespace("dbscan"))
    {
      message('<< Installing dbscan >>')
      install.packages('dbscan', repos = "https://cloud.r-project.org")
    }
    if(!requireNamespace("compaRe"))
    {
      devtools::install_github(repo = 'morchalabi/compaRe', ref = 'master', dependencies = T, force = T)
    }

Report any bug/issue on the [Issue page](https://github.com/morchalabi/compaRe/issues).
