writeSBML <- function(modelData, modelID = "model", outputFile, boundary = "b") {
  if(class(modelData) == "data.frame"){
    # Check valid structure, column names and valid ID's
    modelData <- validateData(modelData = modelData)
    # Remove comments
    modelData <- removeComments(modelData = modelData)
    # Validate stoichiometric syntax
    modelData <- modelData[validateSyntax(modelData[["REACTION"]]),]
    # Convert model to a list
    modelData <- extractData(inputData = modelData)  
  } else {
    stop("Input format not supported.")
  }
  
  # Write SBML model
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
    c(comp, as.vector(sapply(modelData[["COMPARTMENTS"]], function(compartment) {
      paste0('\t\t\t<compartment id="',
             compartment,
             '" name="',
             compartment,
             '"/>')
    })))
  comp <- c(comp, '\t\t</listOfCompartments>')
  mets <- '\t\t<listOfSpecies>'
  mets <-
    c(mets, sapply(modelData[["METABOLITES"]], function(metabolite) {
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
    c(react, unlist(sapply(seq_along(modelData[["REACTIONS"]]), function(reaction) {
      c(
        paste0(
          '\t\t\t<reaction id="R_',
          modelData[["REACTIONS"]][[reaction]][["id"]],
          '"  reversible="',
          modelData[["REACTIONS"]][[reaction]][["reversible"]],
          '">'
        ),
        paste0('\t\t\t\t<notes>'),
        paste0(
          '\t\t\t\t\t<html xmlns="http://www.w3.org/1999/xhtml">',
          if (modelData[["REACTIONS"]][[reaction]][["gpr"]] != "") {
            paste0('<p>GENE_ASSOCIATION: ',
                   modelData[["REACTIONS"]][[reaction]][["gpr"]],
                   '</p>')
          }
          ,'</html>'),
        paste0('\t\t\t\t</notes>'),
        paste0('\t\t\t\t<listOfReactants>'),
        sapply(modelData[["REACTIONS"]][[reaction]][["reactants"]], function(metabolite) {
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
        sapply(modelData[["REACTIONS"]][[reaction]][["products"]], function(metabolite) {
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
          modelData[["REACTIONS"]][[reaction]][["lowbnd"]],
          '"/>'
        ),
        paste0(
          '\t\t\t\t\t\t<parameter id="UPPER_BOUND" value="',
          modelData[["REACTIONS"]][[reaction]][["upbnd"]],
          '"/>'
        ),
        paste0(
          '\t\t\t\t\t\t<parameter id="OBJECTIVE_COEFFICIENT" value="',
          modelData[["REACTIONS"]][[reaction]][["objective"]],
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
