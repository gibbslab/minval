minval : an R package for MINimal VALidation of stoichiometric reactions
======
The **MINVAL** package was designed as a tool to identify orphan metabolites and evaluate the mass and charge balance of stoichometric reactions. **MINVAL** also includes functions to write models in TSV and SBML formats, extract all reactants, products, metabolite names and compartments from a metabolic reconstruction. 

Install:
--------
This package required R version 2.10 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

The official release of minval is available on CRAN. To install from CRAN, use the following command:
```
install.packages("minval", dependencies=TRUE)
```
If you have devtools installed, install the latest stable version this package directly from GitHub:

```
# Install 'devtools' R package
install.packages("devtools")

# LINUX users must install 'libxml' before install 'minval'. Just open a terminal and type:
sudo apt-get install libxml2-dev

# Install 'minval' package
devtools::install_github("gibbslab/minval")
library("minval"")
```

Available functions:
-------------------
|Function | Description |
|:--------|:------------|


Citation
--------
Daniel Osorio, Janneth Gonzalez and Andres Pinzon-Velasco (2016). **minval: MINimal VALidation for Stoichiometric Reactions**. R package version 1.0. https://CRAN.R-project.org/package=minval
