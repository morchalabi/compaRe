# compaRe
To measure similarity between samples, we developed COMPARE. It measures the similarity between two datasets with any number of dimensions and observations.
# Installation
if(!require(devtools))
{
  install.packages("devtools") # If not already installed
}
devtools::install_github(repo = 'morchalabi/compaRe', ref = 'dev', dependencies = T)
