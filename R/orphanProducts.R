#' @export orphanProducts
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

orphanProducts <- function(reactionList, byCompartment=FALSE){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[isValidSyntax(reactionList)]
  # Extract all reactants
  reactant <- unique(unlist(reactants(reactionList)))
  # Extract all products
  product <- unique(unlist(products(reactionList)))
  # Possible candidates to be introduced into the system by exchange reactions or by adding more internal reactions.
  orphan <- rowSums(stoichiometricMatrix(reactionList)!=0)
  orphan <- c(names(orphan)[orphan<2],product[!(product%in%reactant)])
  orphan <- unique(orphan[!orphan%in%reactant[!reactant%in%product]])
  if(length(orphan)==0){
    return (NA)
  } else {
    if (byCompartment == TRUE){
      # Return orphans by compartment
      sapply(compartments(orphan), function(compartment){orphan[grep(paste0("\\[",compartment,"\\]"),orphan)]}, simplify = FALSE)
    } else {
      # Return all reactants never produced in any reaction.
      return(orphan)
    }
  }
}



