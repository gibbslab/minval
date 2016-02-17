minval
======
A package to evaluate the mass balance and extract all reactants, products, orphan metabolites, metabolite names and compartments of stoichiometric reactions. Also are included some options to check the compound names associated to ChEBI database.

Install:
--------
This package required R version 2.10 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

The official release of minval is available on CRAN. To install from CRAN, use the following command:
```
install.packages("minval", dependencies=TRUE)
```
If you have devtools installed, install the latest stable version this package directly from GitHub:

```
library(devtools)
install_github("minval", "dosorio")
library(minval)
```
Available functions:
-------------------
|Function | Description |
|:--------|:------------|
|chebi.candidates|Returns the possible ChEBI names based on compound synonyms|
|chebi.formula|Returns the molecular formula associated to a ChEBI compound name|
|chebi.id|Returns the ChEBI id asociated to a compound name|
|compartments|Identifies the compartments for a set of metabolites|
|is.chebi|Evaluates if a compound name is a ChEBI name|
|is.validsyntax|Evaluates if a stoichiometric reaction has a valid syntax|
|metabolites|Identifies the list of unique metabolites for a set of stoichiometric reactions|
|orphan.products|Identifies the orphan products for a set of stoichometric reactions|
|orphan.reactant|Identifies the orphan reactants for a set of stoichometric reactions|
|products|Identifies the products for a stoichometric reaction|
|reactants|Identifies the reactants for a stoichometric reaction|
|toChEBI|Translates compounds names to ChEBI ids or molecular formulas in a stoichiometric reaction|
|unbalances|Evaluates if a stoichiometric reaction is mass-balanced|

Available datasets
-------------------
| Code        | Description |
|:----------- |:------------|
|chebi|Public data from ChEBI database release 136|
|glugln|Stoichiometric reactions from the reconstructon of the glutamate/glutamine cycle|

Citation
--------
Daniel Osorio, Janneth Gonzalez and Andres Pinzon-Velasco (2016). **minval: MINimal VALidation for Stoichiometric Reactions**. R package version 0.1. https://CRAN.R-project.org/package=minval

