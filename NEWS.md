NEWS
====

#### **minval v.0.8-1**

* An error with the downloadChEBI function when it was used in Windows OS was solved.
* The citation file was added.
* The link with sybilSBML was removed due the package was removed from CRAN

#### **minval v.0.8**

* validateSyntax function now returns FALSE values instead error return
* Minor changes in writeSBMLmod function were made

#### **minval v.0.7**

* A bug related to unmapped metabolites in checkBalance function was solved

#### **minval v.0.6**

* downloadChEBI function was added. This function can download any user-defined release of the ChEBI database.
* checkBalance function now includes support for decimal coefficients in stoichiometric reactions.
* characterizeReactions function was added.
* The RECON 2.04 and an unbalanced Glycolysis model in a human-readable format were included.
* The writeTSV and writeSBML functions were rewritten, the XML package dependency was removed and the support to convert modelOrg objects was added.
* Two vignettes were included
* The orphanMetabolites function was included and the orphanReactants and orphanProducts functions were added as wrappers.

#### **minval v.0.5**

* All functions were vectorized
* Functions to write TSV and SBML models were added
* A function to download the different releases of the ChEBI database was added

#### **minval v.0.4**

* A bug in reactants function was solved. It also affected the orphan.reactants function.

#### **minval v.0.3**
* A bug in products function was solved
* A bug related with spaces between metabolites was solved in several functions

#### **minval v.0.2**
* Some functions were rewrited
* A bug in the compartments function was solved.

#### **minval v.0.1**

* ChEBI database release 136 was added.
* The stoichiometric reactions from the reconstructon of the glutamate/glutamine cycle were added.
