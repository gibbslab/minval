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
splitFormula <- function(coefficient, chemicalFormula) {
  splitAtoms <-regmatches(chemicalFormula,gregexpr("([A-Z]{1}[a-z]?)([0-9]*)", chemicalFormula))
  atoms <- lapply(splitAtoms,function(splitAtoms){
    sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', splitAtoms)
  })
  atomsNumber <-
    lapply(splitAtoms,function(splitAtoms){
      number <- as.numeric(regmatches(x = splitAtoms, m = gregexpr('[0-9]+', splitAtoms)))
      number[is.na(number)] <- 1
      return(number)
    })
  atomsNumber <- sapply(seq_along(atomsNumber),function(molecule){atomsNumber[[molecule]]*coefficient[[molecule]]})
  tapply(unlist(atomsNumber), unlist(atoms), sum)
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
    if (all(is.na(left[[reaction]])) && all(is.na(right[[reaction]]))) {
      return("No compartmentalized reaction")
    } else if (all(is.na(right[[reaction]]))) {
      return("Exchange reaction")
    } else if (length(left[[reaction]]) > 1 ||
               length(right[[reaction]]) > 1) {
      return("Transport reaction")
    } else if (all.equal(target = left[[reaction]], current = right[[reaction]]) == TRUE) {
      return("Compartmentalized reaction")
    } else {
      return("Transport reaction")
    }
  })
}

# Validate modelData
validateData <- function(modelData) {
  names <- colnames(modelData)
  if (length(grep("^ID$", names, ignore.case = TRUE)) == 0) {
    stop("Reaction ID's not found")
  }
  if (!identical(modelData[, "ID"], unique(modelData[, "ID"]))) {
    stop("Reaction ID's must be unique")
  }
  if (length(grep("^REACTION$", names, ignore.case = TRUE)) == 0) {
    stop("REACTIONS not found")
  }
  if (length(grep("^GPR$", names, ignore.case = TRUE)) == 0) {
    stop("GPR not found")
  }
  if (length(grep("^LOWER.BOUND$", names, ignore.case = TRUE)) == 0) {
    stop("LB not found")
  }
  if (length(grep("^UPPER.BOUND$", names, ignore.case = TRUE)) == 0) {
    stop("UB not found")
  }
  if (length(grep("^OBJECTIVE$", names, ignore.case = TRUE)) == 0) {
    stop("OBJECTIVE not found")
  }
  modelData <- as.data.frame.array(modelData)
  modelData[is.na(modelData)] <- ""
  return(modelData)
}

# Remove comments
removeComments <- function(modelData) {
  comments <- grepl("^#", modelData[, 1])
  if (any(comments)) {
    modelData <- modelData[!comments,]
  }
  return (modelData)
}

# Extract Data
extractData <- function(inputData, boundary = "b") {
  exchange <-
    reactionType(inputData[["REACTION"]]) == "Exchange reaction"
  if (any(exchange) == TRUE) {
    inputData[["REACTION"]][exchange] <-
      as.vector(sapply(metabolites(inputData[["REACTION"]][exchange]), function(metabolite) {
        paste0(metabolite,
               " <=> ",
               paste0(
                 metabolites(reactionList = metabolite, woCompartment = TRUE),
                 "[",
                 boundary,
                 "]"
               ))
      }))
  }
  data <- list()
  data$COMPARTMENTS <- compartments(inputData[["REACTION"]])
  data$METABOLITES <-
    metabolites(inputData[["REACTION"]], uniques = TRUE)
  data$REACTIONS <-
    lapply(seq_along(inputData[["REACTION"]]), function(reaction) {
      list(
        id = as.vector(inputData[["ID"]])[reaction],
        reversible = ifelse(
          test = grepl(pattern = "<=>", x = inputData[["REACTION"]][reaction]),
          yes = "true",
          no = "false"
        ),
        gpr = as.vector(inputData[["GPR"]])[reaction],
        reactants = unlist(getLeft(inputData[["REACTION"]][reaction])),
        products = unlist(getRight(inputData[["REACTION"]][reaction])),
        lowbnd = ifelse(
          test = is.numeric(inputData[["LOWER.BOUND"]][reaction]),
          yes = inputData[["LOWER.BOUND"]][reaction],
          no = -1000
        ),
        upbnd = ifelse(
          test = inputData[["UPPER.BOUND"]][reaction] != "",
          yes = inputData[["UPPER.BOUND"]][reaction],
          no = 1000
        ),
        objective = ifelse(
          test = inputData[["OBJECTIVE"]][reaction] != "",
          yes = inputData[["OBJECTIVE"]][reaction],
          no = 0
        )
      )
    })
  return(data)
}

# rearmReactions Rearm the reactions from a stoichiometric matrix
rearmReactions <-
  function(S,
           reversible ,
           type = "SBML",
           boundary = "b") {
    if (type == "SBML") {
      unlist(lapply(seq_len(dim(S)[2]), function(reaction) {
        met = S[, reaction] < 0
        reactants = paste0(abs(S[which(met == TRUE), reaction]), " ", rownames(S)[which(met ==
                                                                                          TRUE)], collapse = " + ")
        met = S[, reaction] > 0
        produts = paste0(abs(S[which(met == TRUE), reaction]), " ", rownames(S)[which(met ==
                                                                                        TRUE)], collapse = " + ")
        paste(reactants,
              produts,
              sep = ifelse(
                test = reversible[reaction],
                yes = " <=> ",
                no = " => "
              ))
      }))
    } else if (type == "TSV") {
      unlist(lapply(seq_len(dim(S)[2]), function(reaction) {
        met = S[, reaction] < 0
        reactants = paste0("(",
                           abs(S[which(met == TRUE), reaction]),
                           ") ",
                           gsub(
                             "\\+",
                             "_ChargedP",
                             gsub("[[:space:]]+", "_", rownames(S)[which(met == TRUE)])
                           ),
                           collapse = " + ")
        if (any(S[, reaction] > 0)) {
          met = S[, reaction] > 0
          produts = paste0("(",
                           abs(S[which(met == TRUE), reaction]),
                           ") ",
                           gsub(
                             "\\+",
                             "_ChargedP",
                             gsub("[[:space:]]+", "_", rownames(S)[which(met == TRUE)])
                           ),
                           collapse = " + ")
        } else {
          met = S[, reaction] > 0
          produts = paste0(abs(S[which(met == TRUE), reaction]), " ", gsub(
            "\\+",
            "_ChargedP",
            gsub("[[:space:]]+", "_", rownames(S)[which(met == TRUE)])
          ), collapse = " + ")
        }
        
        paste(reactants,
              produts,
              sep = ifelse(
                test = reversible[reaction],
                yes = " <==> ",
                no = " --> "
              ))
      }))
    }
  }

# ConvertData from modelOrg for the write functions
convertData <- function(model) {
  data <- NULL
  data$ID <- model@react_id
  data$DESCRIPTION <- model@react_name
  S <- as.matrix(model@S)
  rownames(S) <- model@met_id
  data$REACTION <-
    rearmReactions(S = S, reversible = model@react_rev)
  data$GPR <- model@gpr
  data$LOWER.BOUND <- model@lowbnd
  data$UPPER.BOUND <- model@uppbnd
  data$OBJECTIVE <- model@obj_coef
  data <- data.frame(data)
  return(as.data.frame.array(data))
}
