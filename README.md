# COMPARE-suite
COMPARE-suite is an R package with two functions for similarity measurement and clustering. It employs a dynamic programming algorithm to mesure spatial similarity between two given datasets. One can then use this function to generate a similiarity/affinity matrix which can then be input to the clustering function. The clustering function reports maximal cliques as clusters and dense areas as communities. COMPARE
1. is data-size independent
1. can effectively circumvent the batch effect
1. does not need  dimension reduction
1. does not need subsampling

Read the integrated help page of the R package for more details about setup and running.

For more information about its algorithms check out COMPARE-suite [similarity](https://github.com/morchalabi/COMPARE-suite/wiki/COMPARE-Suite#similarity-matrix-generator) and [clustering](https://github.com/morchalabi/COMPARE-suite/wiki/COMPARE-Suite#clustering) modules. COMPARE-suite is a pipeline for high-throughput screening built on COMPARE.
