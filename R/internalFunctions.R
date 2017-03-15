# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Remove initial and final spaces for a metabolite
removeSpaces <- function(metabolite) {
  metabolite <- gsub(pattern = "^[[:space:]]+",
                     replacement = "",
                     x = metabolite)
  metabolite <- gsub(pattern = "[[:space:]]+$",
                     replacement = "",
                     x = metabolite)
  return(metabolite)
}

# Remove coefficient for a metabolite
removeCoefficients <- function(metabolite) {
  metabolite <-
    gsub(pattern = "^[[:digit:]]+[[:punct:]]*[[:digit:]]*[[:blank:]]+",
         replacement = "",
         x = metabolite)
  return(metabolite)
}

# Remove compartments for a metabolite
removeCompartment <- function(metabolite, rmCoef = FALSE) {
  metabolite <- removeSpaces(metabolite = metabolite)
  if (rmCoef == TRUE) {
    metabolite <- removeCoefficients(metabolite)
  }
  metabolite <-
    gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$", "", metabolite)
  return(metabolite)
}

# Extract metabolites coefficient
coefficients <- function(metabolite) {
  metabolite <- regmatches(
    x = metabolite,
    m = gregexpr(pattern = '^[[:digit:]][[:punct:]]*[[:digit:]]*[[:blank:]]+',
                 text = metabolite)
  )
  metabolite[lengths(metabolite) == 0] <- 1
  metabolite <- as.numeric(removeSpaces(metabolite))
  return(metabolite)
}

# Extract metabolites of the left side in a  biochemical reaction
getLeft <- function(reaction) {
  sides <- unlist(strsplit(x = reaction,split = "[[:space:]]*<?=>[[:space:]]*"))
  metabolites <- removeSpaces(unlist(strsplit(x = sides[1],split = "[[:space:]]+\\+[[:space:]]+")))
  return(metabolites)
}

# Extract metabolites of the right side in a  biochemical reaction
getRight <- function(reaction) {
  sides <- unlist(strsplit(x = reaction,split = "[[:space:]]*<?=>[[:space:]]*"))
  metabolites <- removeSpaces(unlist(strsplit(x = sides[2],split = "[[:space:]]+\\+[[:space:]]+")))
  return(metabolites)
}

