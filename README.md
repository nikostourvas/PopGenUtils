# *PopGenUtils*
  
The *PopGenUtils* package is developed as a way to document a collection of 
helper functions for the analysis of population genetics data.

What can this package do?

- It can gather exported genotyping data from GeneMapper software and generate 
input files for downstream analysis software such as GenAlEx or adegenet.
- It can calculate Probability of Identity (PID) & PIDsibs, which are 
useful statistics in forensic DNA analysis and in cultivar identification (see https://github.com/nikostourvas/PopGenUtils/issues/2#issuecomment-2833505873).
- It provides a wrapper function for the calculation of Polymorphic Information 
Content via the package polysat for genind objects.
- It provides an easy way to create a TreeMix input file for microsatellite data.

You can find more info on how to use the package by typing in the R console:

```
library(PopGenUtils)
help(package = "PopGenUtils")
```

## To install the development version of *PopGenUtils*

You will need the package *devtools*  to be able to install *PopGenUtils*. To install *devtools*:

```
install.packages("devtools")
```

To install *PopGenUtils*:

```
library(devtools)
install_github("nikostourvas/PopGenUtils")
library("PopGenUtils")
```

## Citation
Tourvas N, Ganopoulos I, Koubouris G, et al (2023) Wild and cultivated olive tree genetic diversity in Greece: a diverse resource in danger of erosion. Front Genet 14:. https://doi.org/10.3389/fgene.2023.1298565
