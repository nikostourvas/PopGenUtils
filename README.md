# *PopGenUtils*

The *PopGenUtils* package is developed as a way to document a collection of helper functions for the analysis of population genetics data. Currently at ALPHA stage.

What can this package do?

- It can be used to create a GenAlEx formatted file or a genind object from multiple GeneMapper txt files  
- It can also calculate Probability of Identity (PID) & PIDsibs, which are useful statistics in forensic DNA analysis.


## To install the development version of *PopGenUtils*

You will need the package *devtools*  to be able to install the devel version of *PopGenUtils*. To install *devtools*:

```
install.packages("devtools")
```

To install *PopGenUtils* devel:

```
library(devtools)
install_github("nikostourvas/PopGenUtils")
library("PopGenUtils")
```
