# isBalanced
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

isBalanced <- function(reactionList, referenceData, ids="NAME",mFormula=NULL,mWeight=NULL,mCharge=NULL){
  reactionList <- reactionList[isValidSyntax(as.vector(reactionList))]
  referenceData <- as.data.frame.array(referenceData)
  metabolites <- metabolites(reactionList,woCompartment = TRUE)
  if (is.null(c(mFormula,mWeight,mCharge))){
    stop("Any of mFormula, mWeight or mCharge must be given")
  }
  if(length(c(mFormula,mWeight,mCharge))>1){
    stop("Just one value (mFormula, mWeight or mCharge) must be given to evaluate the reaction balance")
  }
  if(!is.null(mFormula)){
    referenceData <- as.data.frame.array(unique(referenceData[referenceData[,ids]%in%metabolites,c(ids,mFormula)]))
  }
  if(!is.null(mWeight)){
    referenceData <- as.data.frame.array(unique(referenceData[referenceData[,ids]%in%metabolites,c(ids,mWeight)]))
  }
  if(!is.null(mCharge)){
    referenceData <- as.data.frame.array(unique(referenceData[referenceData[,ids]%in%metabolites,c(ids,mCharge)]))
  }
  rownames(referenceData) <- referenceData[,1]
  reactants <- .get.left(reactionList)
  reactants <- lapply(reactants, function(reactants){rep(metabolites(reactants),.coefficients(reactants))})
  products <- .get.right(reactionList)
  products <- lapply(products, function(products){rep(metabolites(products),.coefficients(products))})
  
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