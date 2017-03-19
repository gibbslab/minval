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
      chemicalFormula,
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

# writeSBML
# Write SBML files from a modelData list
# Daniel Osorio <dcosorioh@unal.edu.co>
writeSBML <- function(modelData, modelID, outputFile, boundary = "b") {
  header <- c(
    '<?xml version=\"1.0\" encoding=\"UTF-8\"?>',
    '<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">',
    paste0(
      '\t<model id="',
      paste(modelID, date(), sep = " - "),
      '" name="',
      modelID,
      '">'
    ),
    '\t\t<notes>',
    '\t\t\t<body xmlns="http://www.w3.org/1999/xhtml">',
    '\t\t\t<p> Generated with MINVAL: an R package for MINimal VALidation of stoichiometric reactions </p>',
    '\t\t\t</body>',
    '\t\t</notes>'
  )
  comp <- '\t\t<listOfCompartments>'
  comp <-
    c(comp, sapply(modelData$compartments, function(compartment) {
      paste0('\t\t\t<compartment id="',
             compartment,
             '" name="',
             compartment,
             '"/>')
    }))
  comp <- c(comp, '\t\t</listOfCompartments>')
  mets <- '\t\t<listOfSpecies>'
  mets <-
    c(mets, sapply(modelData$metabolites, function(metabolite) {
      paste0(
        '\t\t\t<species id="M_',
        metabolite,
        '" name="',
        metabolites(metabolite, woCompartment = TRUE),
        '" compartment="',
        compartments(metabolite),
        '" boundaryCondition="',
        ifelse(
          test = compartments(metabolite) == boundary,
          yes = "true",
          no = "false"
        ),
        '"/>'
      )
    }))
  mets <- c(mets, '\t\t</listOfSpecies>')
  react <- '\t\t<listOfReactions>'
  react <-
    c(react, unlist(sapply(seq_along(modelData$reaction), function(reaction) {
      c(
        paste0(
          '\t\t\t<reaction id="R_',
          modelData$reaction[[reaction]]$id,
          '"  reversible="',
          modelData$reaction[[reaction]]$reversible,
          '">'
        ),
        paste0('\t\t\t\t<notes>'),
        paste0(
          '\t\t\t\t\t<html xmlns="http://www.w3.org/1999/xhtml">',
          if (modelData$reaction[[reaction]]$gpr != "") {
            paste0('<p>GENE_ASSOCIATION: ',
                   modelData$reaction[[reaction]]$gpr,
                   '</p>')
          }
          ,'</html>'),
        paste0('\t\t\t\t</notes>'),
        paste0('\t\t\t\t<listOfReactants>'),
        sapply(modelData$reaction[[reaction]]$reactants, function(metabolite) {
          paste0(
            '\t\t\t\t\t<speciesReference species="M_',
            metabolites(metabolite),
            '" stoichiometry="',
            coefficients(metabolite),
            '"/>'
          )
        }),
        paste0('\t\t\t\t</listOfReactants>'),
        paste0('\t\t\t\t<listOfProducts>'),
        sapply(modelData$reaction[[reaction]]$products, function(metabolite) {
          paste0(
            '\t\t\t\t\t<speciesReference species="M_',
            metabolites(metabolite),
            '" stoichiometry="',
            coefficients(metabolite),
            '"/>'
          )
        }),
        paste0('\t\t\t\t</listOfProducts>'),
        paste0('\t\t\t\t<kineticLaw>'),
        paste0(
          '\t\t\t\t\t<math xmlns="http://www.w3.org/1998/Math/MathML">'
        ),
        paste0('\t\t\t\t\t\t<ci>FLUX_VALUE</ci>'),
        paste0('\t\t\t\t\t</math>'),
        paste0('\t\t\t\t\t<listOfParameters>'),
        paste0(
          '\t\t\t\t\t\t<parameter id="LOWER_BOUND" value="',
          modelData$reaction[[reaction]]$lowbnd,
          '"/>'
        ),
        paste0(
          '\t\t\t\t\t\t<parameter id="UPPER_BOUND" value="',
          modelData$reaction[[reaction]]$upbnd,
          '"/>'
        ),
        paste0(
          '\t\t\t\t\t\t<parameter id="OBJECTIVE_COEFFICIENT" value="',
          modelData$reaction[[reaction]]$objective,
          '"/>'
        ),
        paste0('\t\t\t\t\t\t<parameter id="FLUX_VALUE" value="0"/>'),
        paste0('\t\t\t\t\t</listOfParameters>'),
        paste0('\t\t\t\t</kineticLaw>'),
        paste0('\t\t\t</reaction>')
      )
    })))
  react <- c(react, '\t\t</listOfReactions>')
  end <- c('\t</model>', '</sbml>')
  model <- c(header, comp, mets, react, end)
  writeLines(text = model, con = outputFile, sep = "\n")
}
