# metabolites
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

metabolites <- function(reactionList , woCompartment = FALSE, uniques=TRUE){
  reaction <- strsplit(as.vector(reactionList)," <?=>[[:blank:]]*")
  reaction <- lapply(reaction, function(reaction){strsplit(unlist(reaction)," \\+ ")})
  reaction <- lapply(reaction, function(reaction){.remove.spaces(unlist(reaction))})
  reaction <- lapply(reaction, function(reaction){.remove.coefficients(reaction)})
  metabolites <- unlist(reaction)
  if (woCompartment == TRUE) {
    metabolites <- .remove.compartment(metabolites)
  }
  if (uniques == TRUE) {
    metabolites <- unique(metabolites)
  }
  return(metabolites)
}
