# isMassBalanced
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

isBalanced <- function(reactionList, referenceData, ids="NAME",mFormula=NULL,mWeight=NULL,mCharge=NULL){
  reactionList <- as.vector(reactionList)
  referenceData <- as.data.frame.array(referenceData)
  metabolites <- metabolites(reactionList,woCompartment = TRUE)
  
  if(!is.null(mFormula)){
    referenceData <- as.data.frame.array(referenceData[referenceData[,ids]%in%metabolites,c(ids,mFormula)])
  }
  if(!is.null(mWeight)){
    referenceData <- as.data.frame.array(referenceData[referenceData[,ids]%in%metabolites,c(ids,mWeight)])
  }
  if(!is.null(mCharge)){
    referenceData <- as.data.frame.array(referenceData[referenceData[,ids]%in%metabolites,c(ids,mCharge)])
  }
  reactants <- .get.left(reactionList)
  products <- .get.right(reactionList)
  reactants <- lapply(reactants, function(reactants){referenceData[referenceData[,1]%in%reactants,2]})
  products <- lapply(products, function(products){referenceData[referenceData[,1]%in%products,2]})
  
  if(!is.null(mFormula)){
    reactants <- lapply(reactants, function(reactants){.formula2matrix(reactants)})
    products <- lapply(products, function(products){.formula2matrix(products)})
  } else {
    reactants <- lapply(reactants, function(reactants){sum(reactants)})
    products <- lapply(products, function(products){sum(products)})
  }
  balanced <- sapply(seq_along(reactants),function(reaction){all.equal(reactants[[reaction]],products[[reaction]],tolerance = 0.01)})
  return(balanced)
}