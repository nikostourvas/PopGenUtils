# *PopGenUtils*

The *PopGenUtils* package is developed as a way to document a collection of helper functions for the analysis of population genetics data. Currently at ALPHA stage.

The only useful function at this point is called genemapper2df. You can check it out by typing 

```
?genemapper2df
```

in the R console

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
