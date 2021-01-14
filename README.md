# Multiscale Comparative Connectomics

[![arXiv](https://img.shields.io/badge/arXiv-2011.14990-red.svg?style=flat)](https://arxiv.org/abs/2011.14990)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Reproducible code for the results shown in our manuscript [*Multiscale Comparative Connectomics*](https://arxiv.org/abs/2011.14990) (`MCC`).

## Reproducing Figures

Scripts to reproduce the figures in `MCC` are organized below.

### Main Results

#### Figure 1: _Overview of the statistical framework for multiscale comparative connectomics._

#### Figure 2: _Average connectomes for each mouse strain with hierarchical structure labels._

#### Figure 3: _Vertex embeddings of the corpus callosum obtained by the omnibus embedding._

#### Figure 4: _Tractograms of the strongest signal edge, vertex, and community._

#### Figure 5: _Whole-brain embeddings of all mouse connectomes in a low-dimensional space._

#### Figure 6: _Our methods uncover more information about network topology than brain anatomy._

### Supplementary Results

#### S. Figure 1: _DCorr is the best univariate k-sample test for detecting signal edges._

#### S. Figure 2: _The omnibus embedding recovers the most meaningful vertex latent positions._

#### S. Figure 3: _Comparison of representations of community topology._

#### S. Figure 4: _Vertex embeddings of weak signal vertices._

#### S. Figure 5: _Detection of the strongest signal vertices._

#### S. Figure 6: _Tractograms of the weakest signal edge, vertex, and community._

## System Requirements

### Python Dependencies

`MCC` depends on Python 3.6+ and the following libraries
```
graspologic > 0.3
hyppo
```

### R Dependencies
`MCC` depends on R <ADD VERSION NUMBER> and the following libraries
```
ComplexHeatmap
<ADD R DEPENDENCIES>
```
