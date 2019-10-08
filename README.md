## Network analysis of promoter interactions reveals the hierarchical differences in genome organisation between human pluripotent states

___


**This repository contains scripts to reproduce data and figures in Chovanec P. & Collier AJ. *et al*.**

## Installing dependencies

This repository is an RStudio project that uses Packrat for package management and reproducibility. The versions included in Packrat are not the same used in the paper. However, initial testing indicates that the results produced with newer package versions are the same.

When you open the project `CHiC_network.Rproj` in RStudio restore the Packrat snapshot:

`packrat::restore()`


## Getting data

Data that does along with the provided scripts can be found on Open Science Framework: https://osf.io/jp29m/  

Place the scripts inside the `PCHi-C analysis` directory.


Python dependencies not installed by default:

`conda install -c bioconda bx-python`  
`conda install -c conda-forge tqdm`


## Making the interaction network

The `1_make_network.R` script creates the promoter interaction network and outputs it in the graphml format. The visulisation of the network is not performed in R. Instead, we use [Gephi](https://gephi.org/) for all the network visulisations.

If importing your own network into Gephi for the first time, it will not have any layout coordinates.

Depending on the size of the network you should adjust `Scaling` and `Gravity` until a desired layout is achieved, after which ForceAtlas2 can be stopped (it does not stop on its own). In our network the edge weight is set as the linear genomic interaction distance, therefore we have to set `Edge Weight Influence` to `0.0`:

<img src=images/gephi_layout.jpg width="350" height="505" />
