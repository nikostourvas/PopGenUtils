# *PopGenUtils*
  
The *PopGenUtils* package is developed as a way to document a collection of 
helper functions for the analysis of population genetics data.

What can this package do?

- It can gather exported genotyping data from GeneMapper software and generate 
input files for downstream analysis software such as GenAlEx or adegenet.
- It can calculate Probability of Identity (PID) & PIDsibs, which are 
useful statistics in forensic DNA analysis and in cultivar identification.
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
