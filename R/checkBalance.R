#' @export checkBalance
#' @importFrom stats formula
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Evaluate the mass or charge balance for a set of stoichiometric reactions
#' @description For a given set of stoichiometric reactions, this function evaluate the mass or charge balance based in a reference data. Return a boolean value \code{'TRUE'} if reaction is balanced. One of \code{'mFormula'}, \code{'mWeight'} or \code{'mCharge'} arguments must be given.
#' @param reactionList A set of stoichiometric reaction with the following characteristics: \itemize{
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
#' @param referenceData A chemical table containing data to evaluate the balance
#' @param ids A mandatory ID of metabolite names column in the referenceData
#' @param mFormula An optional ID of molecular formula column in the referenceData
#' @param mWeight An optional ID of molecular weight column in the referenceData
#' @param mCharge An optional ID of net charge column in the referenceData
#' @return This function returns a boolean value \code{'TRUE'} if reaction is balanced.
#' @examples
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#'
#' # Loading extrernal chemical information
#' chemicalData <- read.csv2(system.file("extdata", "chemData.csv", package = "minval"))
#' head(chemicalData)
#'
#' # Evaluating mass balance
#' checkBalance(
#' reactionList = glycolysis$REACTION,
#' referenceData = chemicalData,
#' ids = "NAME",
#' mFormula = "FORMULA"
#' )

checkBalance <-
  function(reactionList,
           referenceData,
           ids,
           mFormula = NULL,
           mWeight = NULL,
           mCharge = NULL,
           woCompartment = TRUE) {
    # Validate syntax
    reactionList <-
      as.vector(x = reactionList[validateSyntax(reactionList = reactionList)])
    
    # Extract unique metabolites
    metabolites <- metabolites(reactionList, woCompartment = woCompartment)
    
    # Validate options
    if (!any(colnames(referenceData) %in% ids)) {
      stop("Given ids column not found in referenceData")
    }
    if (is.null(c(mFormula, mWeight, mCharge))) {
      stop("Any of 'mFormula', 'mWeight' or 'mCharge' must be given")
    }
    if (length(c(mFormula, mWeight, mCharge)) > 1) {
      stop(
        "Just one value (mFormula, mWeight or mCharge) must be given to evaluate the reaction balance"
      )
    }
    if (!any(colnames(referenceData) %in% c(mFormula, mWeight, mCharge))) {
      stop("Given options were not found in referenceData")
    }
    
    # Filter reference data
    referenceData <- as.data.frame.array(x = referenceData)
    
    # Filter by Formula
    if (isTRUE(match(mFormula, colnames(referenceData)) > 0)) {
      referenceData <-
        as.data.frame.array(unique(referenceData[referenceData[, ids] %in% metabolites, c(ids, mFormula)]))
    }
    
    # Filter by Mass
    if (isTRUE(match(mWeight, colnames(referenceData)) > 0)) {
      referenceData <-
        as.data.frame.array(unique(referenceData[referenceData[, ids] %in% metabolites, c(ids, mWeight)]))
    }
    
    # Filter by Charge
    if (isTRUE(match(mCharge, colnames(referenceData)) > 0)) {
      referenceData <-
        as.data.frame.array(unique(referenceData[referenceData[, ids] %in% metabolites, c(ids, mCharge)]))
    }
    
    # Change row names
    rownames(referenceData) <- referenceData[, 1]
    
    # Split reaction
    reactants <- getLeft(reactionList)
    products <- getRight(reactionList)
    
    # If molecular formula is given, split formulas and add atoms by type
    if (!is.null(mFormula)) {
      reactants <-
        lapply(reactants, function(reactants) {
          splitFormula(coefficients(reactants), referenceData[metabolites(reactionList = reactants, woCompartment = woCompartment),2])
        })
      products <-
        lapply(products, function(products) {
          splitFormula(coefficients(products), referenceData[metabolites(reactionList = products, woCompartment = woCompartment),2])
        })
      balanced <-
        sapply(seq_along(reactants), function(reaction) {
          isTRUE(all.equal(reactants[[reaction]], products[[reaction]]))
        })
    } else {
      # If other type is given, add the associated values
      reactants <-
        lapply(reactants, function(reactants) {
          sum((referenceData[metabolites(reactionList = reactants, woCompartment = woCompartment),2]) * coefficients(reactants))
        })
      products <-
        lapply(products, function(products) {
          sum((referenceData[metabolites(reactionList = products, woCompartment = woCompartment),2]) * coefficients(products))
        })
      # Check balance
      balanced <-
        sapply(seq_along(reactionList), function(reaction) {
          isTRUE(all.equal(reactants[[reaction]], products[[reaction]], tolerance = 10 ^
                             -min(nchar(
                               sub("[[:digit:]]+\\.", "", referenceData[, 2])
                             ))))
        })
    }
    # Exceptions
    if (dim(referenceData)[1] != length(metabolites)) {
      message(
        "Not all metabolites were mapped correctly, involved reactions will be returned as FALSE"
      )
      balanced[grep(metabolites[!metabolites %in% referenceData[, 1]], reactionList, fixed = TRUE)] <-
        FALSE
    }
    rType <- (reactionType(reactionList) == "Exchange reaction")
    if (any(rType)) {
      message("Exchange reactions identified. It will be returned as TRUE")
      balanced[rType] <- TRUE
    }
    return(balanced)
  }
