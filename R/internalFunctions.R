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
getLeft <- function(reactionList) {
  reactionList <- as.vector(reactionList)
  sides <-
    strsplit(x = reactionList, split = "[[:space:]]*<?=>[[:space:]]*")
  metabolites <-
    lapply(sides, function(sides) {
      removeSpaces(unlist(strsplit(x = sides[1], split = "[[:space:]]+\\+[[:space:]]+")))
    })
  return(metabolites)
}

# Extract metabolites of the right side in a  biochemical reaction
getRight <- function(reactionList) {
  reactionList <- as.vector(reactionList)
  sides <-
    strsplit(x = reactionList, split = "[[:space:]]*<?=>[[:space:]]*")
  metabolites <-
    lapply(sides, function(sides) {
      removeSpaces(unlist(strsplit(x = sides[2], split = "[[:space:]]+\\+[[:space:]]+")))
    })
  return(metabolites)
}

# Split chemical formula
splitFormula <- function(chemicalFormula) {
  splitAtoms <-
    unlist(regmatches(
      formula,
      gregexpr("([A-Z]{1}[a-z]?)([0-9]*)", chemicalFormula)
    ))
  atoms <- sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', splitAtoms)
  atomsNumber <-
    as.numeric(regmatches(x = splitAtoms, m = gregexpr('[0-9]+', splitAtoms)))
  atomsNumber[is.na(atomsNumber)] <- 1
  tapply(atomsNumber, atoms, sum)
}

# Identify the reaction type
reactionType <- function(reactionList) {
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Validate syntax
  reactionList <-
    reactionList[validateSyntax(reactionList = reactionList)]
  # Split reaction
  left <- getLeft(reactionList)
  right <- getRight(reactionList)
  left <-
    lapply(left, function(left) {
      compartments(reactionList = unlist(left))
    })
  right <-
    lapply(right, function(right) {
      compartments(reactionList = unlist(right))
    })
  # Define type of reaction
  sapply(seq_along(reactionList), function(reaction) {
    if (length(right[[reaction]]) > 0) {
      if (length(left[[reaction]]) > 1 | length(right[[reaction]]) > 1) {
        return("Transport reaction")
      } else if (all.equal(target = left[[reaction]], current = right[[reaction]]) == TRUE) {
        return("Compartmentalized reaction")
      } else {
        return("Transport reaction")
      }
    } else {
      return ("Exchange reaction")
    }
  })
}
