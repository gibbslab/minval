#' @export writeTSVmod
#' @importFrom utils write.table
#' @author Daniel Camilo Osorio <dcosorioh@tamu.edu>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Write a model in TSV format for the 'sybil' R package
#' @description This function converts a data.frame or a modelOrg model to TSV format for the 'sybil' R package. TSV format require three \code{'.TSV'} output files (\code{'_react.tsv'}, \code{'_met.tsv'}, \code{'_desc.tsv'}).
#' @param modelData A modelOrg or a data.frame object. If a data.frame is given, it must contain following mandatory colnames: \itemize{
#' \item \code{'ID':} A list of single character strings containing the reaction abbreviations, Entries in the field abbreviation are used as reaction ids, so they must be unique.
#' \item \code{'REACTION':} A set of stoichiometric reaction with the following characteristics: \itemize{
#' \item Arrows symbols must be given in the form \code{'=>'} or \code{'<=>'}
#' \item Inverse arrow symbols \code{'<='} or other types as: \code{'-->'}, \code{'<==>'}, \code{'->'} will not be parsed and will lead to errors.
#' \item Arrow symbols and plus signs (\code{+}) must be surrounded by a space character
#' \item Stoichiometric coefficients must be surrounded by a space character and not by parentheses.
#' \item Each metabolite must have only one stoichiometric coefficient, substituents must be joined to metabolite name by a hyphen (\code{-}) symbol.
#' \item Exchange reactions must have only one metabolite before arrow symbol
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
#' @param outputFile A writable path for the three \code{'.TSV'} output files.
#' @param boundary A single character string specifying the compartment to be used as boundary
#' @return A set of three \code{'.TSV'} files in a valid format to the 'sybil' R package.
#' @examples 
#' #' # Loading a metabolic model
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#' 
#' \dontrun{
#' # Writing a model in TSV format
#' writeTSVmod(modelData = glycolysis,modelID = "Glycolysis",outputFile = "glycolysis")
#' 
#' 
#' # Writing a modelOrg object in a SBML format
#' ## Loading the sybil R package
#' library(sybil)
#' 
#' ## Loading the data
#' data("Ec_core")
#' 
#' ## Writing the modelOrg object in a SBML format
#' writeTSVmod(modelData = Ec_core,modelID = "E.coli",outputFile = "eColi")
#' }
writeTSVmod <-
  function (modelData,
            modelID = "model",
            outputFile,
            boundary = "b") {
    if (is(modelData, "data.frame")) {
      # Check valid structure, column names and valid ID's
      modelData <- validateData(modelData = modelData)
      # Remove comments
      modelData <- removeComments(modelData = modelData)
      # Validate stoichiometric syntax
      modelData <- modelData[validateSyntax(modelData[["REACTION"]]), ]
    } else if (is(modelData, "modelorg")) {
      modelData <- convertData(model = modelData)
    } else {
      stop("Input format not supported.")
    }
    # writeTSV _react.tsv
    outputData <- NULL
    outputData$abbreviation <- modelData[["ID"]]
    outputData$name <- modelData[["DESCRIPTION"]]
    outputData$equation <-
      rearmReactions(
        S = stoichiometricMatrix(modelData[["REACTION"]]),
        reversible = grepl("<=>", modelData[["REACTION"]]),
        type = "TSV",
        boundary = boundary
      )
    outputData$reversible <-
      ifelse(test = grepl("<=>", modelData[["REACTION"]]),
             yes = "reversible",
             no = "irreversible")
    outputData$compartment <-
      unlist(lapply(modelData[["REACTION"]], function(reaction) {
        paste0(compartments(reaction), collapse = " , ")
      }))
    outputData$lowbnd <- as.character(modelData[["LOWER.BOUND"]])
    outputData$uppbnd <- as.character(modelData[["UPPER.BOUND"]])
    outputData$obj_coef <- as.character(modelData[["OBJECTIVE"]])
    outputData$rule <- modelData[["GPR"]]
    outputData$subsystem <- rep("", length(modelData[["REACTION"]]))
    outputData <- data.frame(outputData)
    write.table(
      outputData,
      file = paste0(outputFile, "_react.tsv"),
      quote = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
    
    # writeTSV _met.tsv
    outputData <- NULL
    outputData$abbreviation <-
      gsub("\\+", "_ChargedP", gsub(
        "[[:space:]]+",
        "_",
        metabolites(modelData[["REACTION"]], woCompartment = TRUE)
      ))
    outputData$name <-
      metabolites(modelData[["REACTION"]], woCompartment = TRUE)
    outputData$compartment <-
      compartments(metabolites(modelData[["REACTION"]]))
    outputData <- data.frame(outputData)
    write.table(
      outputData,
      file = paste0(outputFile, "_met.tsv"),
      quote = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
    
    # writeTSV  _desc.tsv
    outputData <- NULL
    outputData$name <- modelID
    outputData$id <- paste(modelID, date(), sep = " - ")
    write.table(
      outputData,
      file = paste0(outputFile, "_desc.tsv"),
      quote = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }
