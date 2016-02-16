# orphan.reactants
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

orphan.reactants <- function(reactionList, byCompartment=FALSE){
  # Extract all reactants
  reactant <- unique(unlist(sapply(reactionList,reactants)))
  # Extract all products
  product <- unique(unlist(sapply(reactionList,products)))

  if (byCompartment == TRUE){
    # Identifies all reactants never produced in any reaction.
    orphan <- reactant[!reactant%in%product]
    # Return orphans by compartment
    sapply(compartments(orphan), function(comp){orphan[grep(comp,orphan)]}, simplify = FALSE)
  } else {
    # Return all reactants never produced in any reaction.
    # Possible candidates to be introduced into the system by exchange reactions or by adding more internal reactions.
    reactant[!reactant%in%product]
  }
}

