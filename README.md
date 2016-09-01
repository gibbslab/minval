minval : MINimal VALidation of stoichiometric reactions
======
The **MINVAL** package was designed as a tool to identify orphan metabolites and evaluate the mass and charge balance of stoichometric reactions. **MINVAL** also includes functions to write TSV and SBML models, extract all reactants, products, metabolite names and compartments from a metabolic reconstruction. 

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
|compartments|Extract the list of unique compartments for the metabolites of a set of stoichiometric reactions.|
|convert2sbml|Write a SBML file.|
|convert2sbmlR|Convert a data.frame data to a SBMLR object|
|convert2tsv|Write a TSV model for the 'sybil' R package.|
|getChEBI|Download the ChEBI database|
|isBalanced||
|isValidSyntax||
|metabolites||
|orphanProducts||
|orphanReactants||
|products||
|reactants||
|stoichiometricMatrix||
|xls2sbml||
|xls2tsv||

Citation
--------
Daniel Osorio, Janneth Gonzalez and Andres Pinzon-Velasco (2016). **minval: MINimal VALidation for Stoichiometric Reactions**. R package version 0.5. https://CRAN.R-project.org/package=minval
