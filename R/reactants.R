# reactants | Identifies the reactants for a stoichometric reaction
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

reactants <- function(reaction){
  reaction <- unlist(lapply(strsplit(as.vector(reaction)," => "),function(reaction){reaction[[1]]}))
  reaction <- strsplit(reaction, " <=> ")
  # Split independient reactants
  reaction <- lapply(reaction, function(reaction){strsplit(reaction," \\+ ")})
  # Remove spaces and report uniques
  reaction <- lapply(reaction, function(reaction){unique(.remove.spaces(unlist(reaction)))})
  # Use a regex to extract stoichiometric coefficients and separate the metabolite name
  reactants <- lapply(reaction, function(reaction){unique(.remove.coefficients(reaction))})
  if (length(reactants)==1){
    return(unlist(reactants))
  } else {
    return(reactants)
  }
}
