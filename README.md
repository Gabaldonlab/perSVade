# perSVade: personalized Structural Variation detection
![Latest Version](https://img.shields.io/github/v/tag/gabaldonlab/perSVade?label=Latest%20Version)
[![DOI](https://zenodo.org/badge/263600715.svg)](https://zenodo.org/badge/latestdoi/263600715)
![Docker Pulls](https://img.shields.io/docker/pulls/mikischikora/persvade)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/repository/docker/mikischikora/persvade)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

perSVade runs structural variation (SV), small variant (SNPs and IN/DELs) and read depth-based Copy Number Variation (CNV) calling and annotation for WGS datasets. The only required input is a set of paired-end short reads and a reference genome. Everything with a few simple commands.

The SV calling pipeline finds breakpoints with  [GRIDSS](https://github.com/PapenfussLab/gridss) and summarizes them into complex structural variants with [CLOVE](https://github.com/PapenfussLab/clove), with some added features. perSVade provides an automated benchmarking and parameter selection for these methods in any genome or sequencing run. This is useful for species without reccommended running and filtering parameters. In addition, it provides an automated report of the SV calling accuracy on simulations, useful to assess the confidence of the results. 

Check the [wiki](https://github.com/Gabaldonlab/perSVade/wiki) in order to install and use the pipeline. Also, check our [video tutorial](https://www.youtube.com/watch?v=Nve8LGJ7AC0) for a quick introduction to perSVade.
