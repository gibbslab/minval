#' @export writeSBML
#' @author  Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Write a model in SBML format
#' @description This function converts a data.frame or a modelOrg object to a valid SBML file. The Systems Biology Markup Language (SBML) is a representation format, based on XML, for communicating and storing computational models of biological processes.
#' @param modelData A modelOrg or a data.frame object. If a data.frame is given, it must contain following mandatory colnames: \itemize{
#' \item \code{'ID':} A list of single character strings containing the reaction abbreviations, Entries in the field abbreviation are used as reaction ids, so they must be unique.
#' \item \code{'REACTION':} A set of stoichiometric reaction with the following characteristics: \itemize{
#' \item Arrows symbols must be given in the form \code{'=>'} or \code{'<=>'}
#' \item Inverse arrow symbols \code{'<='} or other types as: \code{'-->'}, \code{'<==>'}, \code{'->'} will not be parsed and will lead to errors.
#' \item Arrow symbols and plus signs (\code{+}) must be surrounded by a space character
#' \item Stoichiometric coefficients must be surrounded by a space character and not by parentheses.
#' \item Each metabolite must have only one stoichiometric coefficient, substituents must be joined to metabolite name by a hyphen (\code{-}) symbol.
#' \item Exchange reactions have only one metabolite before arrow symbol
#' \item Compartments must be given between square brackets ([compartment]) joined at the end of metabolite name
#' }
#' Some examples of valid stoichiometric reactions are: \itemize{
#' \item \code{H2O[c] + Urea-1-Carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]}
#' \item \code{ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]}
#' \item \code{CO2[c] <=> }
#' }
#' \item \code{'GPR':} A set of genes joined by boolean operators as AND or OR, rules may be nested by parenthesis. (optional: column can be empty),
#' \item \code{'LOWER.BOUND':} A list of numeric values containing the lower bounds of the reaction rates.
#' If not set, zero is used for an irreversible reaction and -1000 for a reversible reaction. (optional: column can be empty),
#' \item \code{'UPPER.BOUND':} A list of numeric values containing the upper bounds of the reaction rates.
#' If not set, 1000 is used by default. (optional: column can be empty),
#' \item \code{'OBJECTIVE':} A list of numeric values containing objective values (-1, 0 or 1) for each reaction (optional: column can be empty).
#' }
#' @param modelID A single character string giving the modelID
#' @param outputFile A writable path for the output 'SBML' file to be generate
#' @param boundary A single character string specifying the compartment to be used as boundary, \code{'b'} is used by default
#' @return A valid SBML file
#' @examples
#' # Loading a metabolic model
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#'
#' # Writing a model in SBML format
#' writeSBML(modelData = glycolysis,modelID = "Glycolysis",outputFile = "glycolysis.xml")
#' 
#' \dontrun{
#' # Writing a modelOrg object in a SBML format
#' ## Loading the sybil R package
#' library(sybil)
#' 
#' ## Loading the data
#' data("Ec_core")
#' 
#' ## Writing the modelOrg object in a SBML format
#' writeSBML(modelData = Ec_core,modelID = "E.coli",outputFile = "eColi.xml")
#' }

writeSBML <-
  function(modelData,
           modelID = "model",
           outputFile,
           boundary = "b") {
    if (class(modelData) == "data.frame") {
      # Check valid structure, column names and valid ID's
      modelData <- validateData(modelData = modelData)
      # Remove comments
      modelData <- removeComments(modelData = modelData)
      # Validate stoichiometric syntax
      modelData <- modelData[validateSyntax(modelData[["REACTION"]]), ]
    } else if (class(modelData) == "modelorg") {
      modelData <- convertData(model = modelData)
    } else {
      stop("Input format not supported.")
    }
    # Convert model to a list
    modelData <- extractData(inputData = modelData)
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
            ,
            '</html>'
          ),
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
    model <- as.vector(c(header, comp, mets, react, end))
    writeLines(text = model, con = outputFile, sep = "\n")
  }
