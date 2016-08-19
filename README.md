minval : MINimal VALidation of stoichiometric reactions
======
The **minval** package was designed as a tool to identify orphan metabolites and the mass unbalanced reactions in a set of stoichometry reactions, it also permits to extract all reactants, products, metabolite names and compartments from a metabolic reconstruction.  Moreover specific functions to map compound names associated to the Chemical Entities of Biological Interest (ChEBI) database are also included.

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

# Install 'minval' package
devtools::install_github("gibbslab/minval")
library(minval)
```

Available functions:
-------------------
|Function | Description |
|:--------|:------------|
|compartments|Identifies the compartments for a set of metabolites|
|convert2SBML|Converts a data.frame to a SBML file
|convert2SBMLR|Converts a data.frame to a SBMLR object
|convert2TSV|Converts a data.frame to a TSV model for 'sybil' package
|is.validSyntax|Evaluates if a stoichiometric reaction has a valid syntax|
|metabolites|Identifies the list of unique metabolites for a set of stoichiometric reactions|
|orphanProducts|Identifies the orphan products for a set of stoichometric reactions|
|orphanReactant|Identifies the orphan reactants for a set of stoichometric reactions|
|products|Identifies the products for a stoichometric reaction|
|reactants|Identifies the reactants for a stoichometric reaction|
|stoichiometricMatrix|Generate the stoichiometric matrix from a set of stoichiometric reactions|

Available datasets
-------------------
| Code        | Description |
|:----------- |:------------|
|chebi|Public data from ChEBI database release 136|
|glugln|Stoichiometric reactions from the reconstructon of the glutamate/glutamine cycle|

Citation
--------
Daniel Osorio, Janneth Gonzalez and Andres Pinzon-Velasco (2016). **minval: MINimal VALidation for Stoichiometric Reactions**. R package version 0.1. https://CRAN.R-project.org/package=minval
