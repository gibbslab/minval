# reactants | Identifies the reactants for a stoichometric reaction
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

reactants <- function(reactionList){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[is.validSyntax(reactionList)]
  # Extract reactants for irreversible reactions
  reaction <- unlist(lapply(strsplit(reactionList,"[[:blank:]]+=>[[:blank:]]+"),function(reactionList){reactionList[[1]]}))
  # Extract metabolites for reversible reactions
  reaction <- strsplit(reaction, "[[:blank:]]+<=>[[:blank:]]+")
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
