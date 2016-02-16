# products
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

products <- function(rxn){
  # Identifies if stoichiometric reaction is not reversible. In this case:
  if (length(grep("=",sub("(.*) <=> (.*)","\\2",rxn)))>0){
    # Use a regex to extract stoichiometric coefficients and separate the metabolite name
    # Return all the products
    sub("^[[:digit:]]?.?[[:digit:]]+ ","\\2",strsplit(sub("(.*) => (.*)","\\2",rxn),fixed = TRUE,split = " + ")[[1]])
  }
  # In contrary case:
  else {
    # Use a regex to extract stoichiometric coefficients and separate the metabolite name
    # Return both products and reactants
    c(sub("^[[:digit:]]?.?[[:digit:]]+ ","\\2",strsplit(sub("(.*) <=> (.*)","\\1",rxn),fixed = TRUE,split = " + ")[[1]]),sub("^[[:digit:]]?.?[[:digit:]]+ ","\\2",strsplit(sub("(.*) <=> (.*)","\\2",rxn),fixed = TRUE,split = " + ")[[1]]))
  }
}

