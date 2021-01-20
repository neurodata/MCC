# Multiscale Comparative Connectomics

[![arXiv](https://img.shields.io/badge/arXiv-2011.14990-red.svg?style=flat)](https://arxiv.org/abs/2011.14990)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Reproducible code for the results shown in our manuscript [*Multiscale Comparative Connectomics*](https://arxiv.org/abs/2011.14990) (`MCC`).

> Vivek Gopalakrishnan, Jaewon Chung, Eric Bridgeford, Benjamin D. Pedigo, Jesús Arroyo, Lucy Upchurch, G. Allan Johnson, Nian Wang, Youngser Park, Carey E. Priebe, and Joshua T. Vogelstein. “Multiscale Comparative Connectomics”. arXiv:2011.14990 (Nov. 2020).

## Table of Figures and Tables

#### [Figure 1: _Overview of the statistical framework for multiscale comparative connectomics._](#figure-1)

#### [Figure 2: _Average connectomes for each mouse strain with hierarchical structure labels._](#figure-2)

#### [Figure 3: _Vertex embeddings of the corpus callosum obtained by the omnibus embedding._](#figure-3)

#### [Figure 4: _Tractograms of the strongest signal edge, vertex, and community._](#figure-4)

#### [Figure 5: _Whole-brain embeddings of all mouse connectomes in a low-dimensional space._](#figure-5)

#### [Figure 6: _Our methods uncover more information about network topology than neuroanatomical metrics._](#figure-6)

## System Requirements

### Python Dependencies

`MCC` depends on Python 3.8 and the following libraries:
```
graspologic @ git+git://github.com/microsoft/graspologic@dev
ipykernels
rpy2
scikit-image
statsmodels
```

### R Dependencies
`MCC` depends on R 3.6.3 and the following libraries:
```
cdcsis
ComplexHeatmap
circlize
data.table
future
future.apply
ggplot2
igraph
mltools
tidyverse
```

## Reproducing Figures

Scripts to reproduce the figures in `MCC` are organized below.

#### Figure 1
![Fig1](code/figures/1_framework.jpg)

Run `code/1_statistical_framework_graphs.ipynb`.
This script uses `igraph` to generate the sample connectomes and graph models seen above.

#### Figure 2
![Fig2](code/figures/2_connectome.jpg)

#### Figure 3

#### Figure 4

#### Figure 5

#### Figure 6
