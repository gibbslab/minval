# metabolites
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

metabolites <- function(reactionList , woCompartment = FALSE){
  # Extract and return the unique metabolites list
  mets <- unique(
    # Remove the list format and convert it to an array
    unlist(
      # Join reactants and products from the resultant lists
      c(
        # Extract all reactants from the stoichiometric reactions
        sapply(reactionList,reactants),
        # Extract all products from the stoichiometric reactions
        sapply(reactionList,products)
      )
    )
  )

  if (woCompartment == TRUE){
    unique(sapply(mets, function(met){.metname(met,rm.coef = TRUE)}))
  } else{
    mets
  }
}
