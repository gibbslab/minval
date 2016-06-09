# products | Identifies the reactants for a stoichometric reaction
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

products <- function(reaction){
  reaction <- strsplit(as.vector(reaction)," => ")
  reaction[lengths(reaction)>1] <- lapply(reaction[lengths(reaction)>1],function(reaction){reaction[[2]]})
  reaction <- strsplit(unlist(reaction), " <=> ")
  # Split independient reactants
  reaction <- lapply(reaction, function(reaction){strsplit(reaction," \\+ ")})
  # Remove spaces and report uniques
  reaction <- lapply(reaction, function(reaction){unique(.remove.spaces(unlist(reaction)))})
  # Use a regex to extract stoichiometric coefficients and separate the metabolite name
  products <- lapply(reaction, function(reaction){unique(.remove.coefficients(reaction))})
  if (length(products)==1){
    return(unlist(products))
  } else {
    return(products)
  }
}

