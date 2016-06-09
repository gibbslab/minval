# metabolites
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

metabolites <- function(reactionList , woCompartment = FALSE){
  reaction <- strsplit(as.vector(reactionList)," <?=> ")
  reaction <- lapply(reaction, function(reaction){strsplit(unlist(reaction)," \\+ ")})
  reaction <- lapply(reaction, function(reaction){unique(.remove.spaces(unlist(reaction)))})
  reaction <- lapply(reaction, function(reaction){unique(.remove_coefficients(reaction))})
  metabolites <- unique(unlist(reaction))
  if (woCompartment == TRUE){
    return(unique(.metname(metabolites)))
  } else{
    return(metabolites)
  }
}
