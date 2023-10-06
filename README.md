# *Generation of de novo miRNAs from template switching during DNA replication*

The code here allows reproducing the computational analyses in the study. Small datasets are provided to replicate specific steps and full data are available elsewhere.



## Instructions for replicating the analyses

Instructions for replicating the analyses are given in the Jupyter Notebook  [doc/project_description.ipynb](doc/project_description.ipynb). Instructions for the installation of the tools and software required for the analyses are given below.
The associated data and results directory can be downloaded from: https://etsin.fairdata.fi/dataset/d3f0360c-cd36-46f5-b5ff-32398281d48d

## Data for TSM-origin microRNAs

Ancestral sequence reconstructions, phylogenetic trees, TSM alignment solutions and predicted secondary structures are provided under [data/microRNAs](data/). Subdirectories are named in format 'EnsembleID_microRnaID', TSM solutions  'EnsembleID_TSM-number' and predicted secondary structures 'EnsembleID_TSM-number_Tree-node-qry/ref/anc'. The sequence containing a TSM is 'qry', and its parent and grand-parent are 'ref' and 'anc', respectively.


## Tools and software

### Conda environment

The required python packages and their dependecies can installed 
by creating a conda environment:

``conda env create -f microRNA_TSM.yml``

The YAML description is given in [scripts/microRNA_TSM.yml](scripts/microRNA_TSM.yml).

For the R code:

``conda env create -f bioconductor.yml ``

The YAML code is given in [scripts/bioconductor.yml](scripts/biocondcutor.yml).

### Software

In addition to the python packages, the following software are required:

- bedtools v.2.26.0 (https://bedtools.readthedocs.io/en/latest/)

- Vienna package v.2.4.14 (https://www.tbi.univie.ac.at/RNA/#)

- MongoDB (https://docs.mongodb.com/guides/server/install/ )

- Dustmasker v.1.0 (https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/)

- Pagan2 (https://github.com/ariloytynoja/pagan2-msa)

- python library for FPA (https://github.com/ariloytynoja/fpa; see below)


The python API for FPA is in the file [scripts/cpp/fpa4.cpp](scripts/cpp/fpa4.cpp). Using GCC, it can be compiled as:

```
c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` fpa4.cpp -o fpa_ext62`python3-config --extension-suffix` 
```

### Tools

Python, bash and R scripts are given in directories [scripts/python/](scripts/python/), [scripts/bash/](scripts/bash/) and [scripts/R/](scripts/R/). 
