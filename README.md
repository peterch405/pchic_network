## Network analysis of promoter interactions reveals the hierarchical differences in genome organisation between human pluripotent states

___


**This repository contains scripts to reproduce the analysis and figures in Chovanec P. & Collier AJ. *et al*.**

Preprint:
**Network analysis of promoter interactions reveals the hierarchical differences in genome organisation between human pluripotent states.** Peter Chovanec, Amanda J. Collier, Christel Krueger, Csilla Várnai, Stefan Schoenfelder, Anne Corcoran, Peter J. Rugg-Gunn
bioRxiv 2019.12.13.875286; doi: https://doi.org/10.1101/2019.12.13.875286 

## Installing dependencies


This repository is an RStudio project that relies on renv for package management to ensure reproducibility.

Ensure that you have renv is installed. It can be installed by running the following command:

`install.packages("renv")`

When you open the project `CHiC_network.Rproj` in RStudio, first restore the renv snapshot by running:

`renv::restore()`

#### Install time

Installation of all the packages takes ~30 minutes on a normal computer.

## Getting data

Data that goes alongside the provided scripts can be found on Open Science Framework: https://osf.io/jp29m/  (The scripts have also been placed on OSF in the correct folders to facilitate re-running)

Place the scripts inside the `PCHi-C analysis` directory.


Python dependencies not installed by default:

`conda install -c bioconda bx-python`  
`conda install -c conda-forge tqdm`


## Making the interaction network

The `1_make_network.R` script creates the promoter interaction network and outputs it in the graphml format. The visulisation of the network is not performed in R. Instead, we use [Gephi](https://gephi.org/) for all the network visulisations.

If importing your own network into Gephi for the first time, it will not have any layout coordinates.

Depending on the size of the network you should adjust `Scaling` and `Gravity` until a desired layout is achieved, after which ForceAtlas2 can be stopped (it does not stop on its own). In our network the edge weight is set as the linear genomic interaction distance, therefore we have to set `Edge Weight Influence` to `0.0`:

<img src=images/gephi_layout.jpg width="350" height="505" />

#### Run time

Running all the scripts to reproduce paper figures takes ~2 hours on a normal computer.
