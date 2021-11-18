# *PopGenUtils*

<!-- badges: start -->
  [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
  <!-- badges: end -->
  
The *PopGenUtils* package is developed as a way to document a collection of 
helper functions for the analysis of population genetics data.

What can this package do?

- It can gather exported genotyping data from GeneMapper software and generate 
input files for downstream analysis software such as GenAlEx or adegenet.
- It can calculate Probability of Identity (PID) & PIDsibs, which are 
useful statistics in forensic DNA analysis.
- It provides a wrapper function for the calculation of Polymorphic Information 
Content via the package polysat for genind objects.

You can find more info on how to use the package by typing in the R console:

```
library(PopGenUtils)
help(package = "PopGenUtils")
```

## To install the development version of *PopGenUtils*

You will need the package *devtools*  to be able to install the devel version 
of *PopGenUtils*. To install *devtools*:

```
install.packages("devtools")
```

To install *PopGenUtils* devel:

```
library(devtools)
install_github("nikostourvas/PopGenUtils")
library("PopGenUtils")
```
