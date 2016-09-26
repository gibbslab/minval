#' @export products
#' @title Identify the reactants of a stoichometric reaction
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

products <- function(reactionList){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Remove reaction with invalid syntax
  reactionList <- reactionList[isValidSyntax(reactionList)]
  # Extract reactants for irreversible reactions
  reaction <- strsplit(reactionList,"[[:blank:]]+=>[[:blank:]]+")
  reaction[lengths(reaction)>1] <- lapply(reaction[lengths(reaction)>1],function(reaction){reaction[[2]]})
  # Extract metabolites for reversible reactions
  reaction <- strsplit(unlist(reaction), "[[:blank:]]+<=>[[:blank:]]+")
  # Split independient reactants
  reaction <- lapply(reaction, function(reaction){strsplit(reaction,"[[:blank:]]+\\+[[:blank:]]+")})
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

