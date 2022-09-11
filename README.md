## Nonparametric estimation of state occupation probabilities from multistate models with current-status data

This repository contains code for computing the state occupation probabilities (SOP) for current-status data from of a general multistate model. The code estimates the  SOP for the setting where the current-status data is either uncorrelated or cluster-correlated. In the case where the current-status data is cluster-correlated, the cluster sizes may be subject to informative cluster size (ICS), the methods implemented in the code also adjusts for the ICS.

In addition, this repository contains the a real-world current-status data set from a study where the periodontal disease health of Gullah-speaking African-Americans diabetics (GAAD) was examined at a given inspection time [[1]](#1). Incidence of periodontal disease was determined by measuring the clinical attachment level at six sites of each tooth, excluding molars, in the subjects' mouth. The dataset comprise both clinical and demographic information for 288 subjects that participated in the study. 

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

## References
<a id="1">[1]</a> 
Fernandes JK, Wiegand RE, Salinas CF, Grossi SG, Sanders JJ, Lopes-Virella MF, Slate EH (2009). 
Periodontal disease status in Gullah African Americans with type 2 diabetes living in South Carolina. 
*Journal of periodontology*, **80**(7), 1062–1068.
