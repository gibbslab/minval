# orphan.products
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

orphan.products <- function(reactionList, byCompartment=FALSE){
  # Extract all reactants
  reactant <- unique(unlist(sapply(reactionList,reactants)))
  # Extract all products
  product <- unique(unlist(sapply(reactionList,products)))

  if (byCompartment == TRUE){
    # Identifies orphan products
    orphan <- product[!product%in%reactant]
    # Return orphans by compartment
    sapply(compartments(orphan), function(comp){orphan[grep(paste("[[:punct:]]",comp,"[[:punct:]]",sep = ""),orphan)]}, simplify = FALSE)
  } else {
    # Return all products that are not consumed in any reaction. Possible dead ends.
    product[!product%in%reactant]
  }
}



