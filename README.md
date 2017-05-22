minval : an R package for MINimal VALidation of stoichiometric reactions
======
The genome-scale metabolic reconstructions, a compilation of all stoichiometric reactions that can describe the entire cellular metabolism of an organism, have become an indispensable tool for our understanding of biological phenomena, covering fields that range from systems biology to bioengineering. Evaluation of metabolic reconstructions are generally carried through Flux Balance Analysis, an optimization method where the biological sense of optimal solution is sensitive to thermodynamic unbalance, caused by the presence of stoichiometric reactions whose compounds are not produced or consumed in any other reaction (orphan metabolites) and by mass unbalanced stoichiometric reactions. The **MINVAL** package was designed as a tool to identify orphan metabolites and evaluate the mass and charge balance of stoichiometric reactions. **MINVAL** package also includes functions to characterize and write models in TSV and SBML formats, extract all reactants, products, metabolite names and compartments from a metabolic reconstruction.

Install:
--------
This package required R version 2.10 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

The official release of **MINVAL** is available on CRAN. To install from CRAN, use the following command:
```{r}
install.packages("minval", dependencies = TRUE)
```
If you have devtools installed, install the latest stable version this package directly from GitHub:

```{r}
# Install 'minval' package
devtools::install_github("gibbslab/minval")
library("minval")
```

Available functions:
-------------------
|Function | Description |
|:--------|:------------|
|characterizeReactions|Characterize stoichiometric reactions by compartments and reaction type|
|checkBalance|Evaluate the mass or charge balance for a set of stoichiometric reactions|
|compartments|Extract the compartments associated to metabolites of a set of stoichiometric reactions|
|downloadChEBI|Download the ChEBI database|
|metabolites|Identify the list of metabolites for a set of stoichiometric reactions|
|orphanMetabolites|Identify the orphan metabolites of a set of stoichiometric reactions|
|products|Identify the products of a stoichometric reaction|
|reactants|Identify the reactants of a stoichometric reaction|
|stoichiometricMatrix|Build the stoichiometric matrix for a set of stoichiometric reactions|
|validateSyntax|Evaluate if a stoichiometric reaction has a valid syntax|
|writeSBMLmod|Write a model in SBML format|
|writeTSVmod|Write a model in TSV format for the 'sybil' R package|

Citation
--------
Daniel Osorio, Janneth Gonzalez and Andres Pinzon (2016). **minval: MINimal VALidation for Stoichiometric Reactions**. R package version 0.6. https://CRAN.R-project.org/package=minval
