## Nonparametric estimation of state occupation probabilities from multistate models with current-status data

This repository contains code for computing the state occupation probabilities (SOP) for current-status data from of a general multistate model. The code estimates the  SOP for the setting where the current-status data is either uncorrelated or cluster-correlated. In the case where the current-status data is cluster-correlated, the cluster sizes may be subject to informative cluster size (ICS), the methods implemented in the code also adjusts for the ICS.

This package uses functions from `graph` and `Rgraphviz` to create graph structure represented in terms of nodes and an edge list for the multistate model.

## Installation
The package depends on `graph` and `Rgraphviz`, these packages can be installed via `BiocManager`:

   if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")

   BiocManager::install("graph")
   BiocManager::install("Rgraphviz")


The `mspack2` package can be installed using the `devtools` packgage:

   library(devtools) 
   devtools::install_github(repo='"samuelanyaso/mspack2") 
