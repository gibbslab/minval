#' @export isBalanced
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Evaluate the mass or charge balance for a set of stoichiometric reactions
#' @description For a given set of stoichiometric reactions, this function evaluate the mass or charge balance based in a reference data. Return a boolean value \code{'TRUE'} if reaction is balanced. One of \code{'mFormula'}, \code{'mWeight'} or \code{'mCharge'} arguments must be given.
#' @param reactionList A set of stoichiometric reaction with the following format: \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param referenceData A chemical table containing data to evaluate the balance
#' @param ids A mandatory id of metabolites id column in the referenceData
#' @param mFormula An optional id of molecular formula column in the referenceData
#' @param mWeight An optional id of molecular weight column in the referenceData
#' @param mCharge An optional id of net charge column in the referenceData
#' @return This function returns a boolean value \code{'TRUE'} if reaction is balanced. Reactions without valid syntax are removed of the reactionList.
#' @keywords check mass charge balance genome scale metabolic reconstruction
#' @examples 
#' # Using external chemical table
#' chemicalData <- read.csv2(system.file("extdata", "chemData.csv", package = "minval"))
#' head(chemicalData)
#' 
#' # Loading stoichiometric reactions
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#' 
#' # Evaluating mass balance
#' isBalanced(
#'  reactionList = glycolysis$REACTION, 
#'  referenceData = chemicalData, 
#'  ids = "NAME", 
#'  mFormula = "FORMULA"
#'  )
#' 
#' # Evaluating charge balance
#' isBalanced(
#'  reactionList = glycolysis$REACTION, 
#'  referenceData = chemicalData, 
#'  ids = "NAME",
#'  mCharge = "CHARGE"
#'  )

isBalanced <- function(reactionList, referenceData, ids, mFormula=NULL, mWeight=NULL, mCharge=NULL){
  reactionList <- as.vector(reactionList[isValidSyntax(as.vector(reactionList))])
  referenceData <- as.data.frame.array(referenceData)
  metabolites <- metabolites(reactionList,woCompartment = TRUE)
  if (is.null(c(mFormula,mWeight,mCharge))){
    stop("Any of mFormula, mWeight or mCharge must be given")
  }
  if(length(c(mFormula,mWeight,mCharge))>1){
    stop("Just one value (mFormula, mWeight or mCharge) must be given to evaluate the reaction balance")
  }
  if(isTRUE(match(mFormula,colnames(referenceData))>0)){
    referenceData <- as.data.frame.array(unique(referenceData[referenceData[,ids]%in%metabolites,c(ids,mFormula)]))
  }
  if(isTRUE(match(mWeight,colnames(referenceData))>0)){
    referenceData <- as.data.frame.array(unique(referenceData[referenceData[,ids]%in%metabolites,c(ids,mWeight)]))
  }
  if(isTRUE(match(mCharge,colnames(referenceData))>0)){
    referenceData <- as.data.frame.array(unique(referenceData[referenceData[,ids]%in%metabolites,c(ids,mCharge)]))
  }
  rownames(referenceData) <- referenceData[,1]
  reactants <- lapply(reactionList, .get.left)
  reactants <- lapply(reactants, function(reactants){rep(metabolites(reactants,woCompartment = TRUE,uniques = FALSE),.coefficients(reactants))})
  products <- lapply(reactionList, .get.right)
  products <- lapply(products, function(products){(rep(metabolites(products,woCompartment = TRUE,uniques = FALSE),.coefficients(products)))})
  
  reactants <- lapply(reactants, function(reactants){referenceData[reactants,2]})
  products <- lapply(products, function(products){referenceData[products,2]})
  if(!is.null(mFormula)){
    reactants <- lapply(reactants, function(reactants){.formula2matrix(reactants)})
    products <- lapply(products, function(products){.formula2matrix(products)})
    balanced <- sapply(seq_along(reactants),function(reaction){isTRUE(all.equal(reactants[[reaction]],products[[reaction]]))})
  } else {
    reactants <- lapply(reactants, function(reactants){sum(as.numeric(reactants))})
    products <- lapply(products, function(products){sum(as.numeric(products))})
    balanced <- sapply(seq_along(reactants),function(reaction){isTRUE(all.equal(reactants[[reaction]],products[[reaction]],tolerance = 10^-min(nchar(sub("[[:digit:]]+\\.","",referenceData[,2])))))})
  }
  if(dim(referenceData)[1]!=length(metabolites)){
    message("Not all metabolites were mapped correctly, involved reactions will be returned as NA")
    balanced[grep(metabolites[!metabolites%in%referenceData[,1]],reactionList,fixed = TRUE)] <- NA
  }
  return(balanced)
}


