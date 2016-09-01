#' @export mapReactions
#' @author Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Return reactions of a reference data by a selected column

mapReactions <- function(reactionList, referenceData, by, inverse=FALSE ){
  if(!is.null(dim(referenceData)) && is.null(dim(reactionList))){
    referenceData <- as.data.frame.array(referenceData)
    if(by%in%colnames(referenceData)){
      if(inverse == FALSE){
        referenceData[referenceData[,match(by,colnames(referenceData))]%in%reactionList,]
      } else {
        referenceData[!referenceData[,match(by,colnames(referenceData))]%in%reactionList,]
      }
    }
  }
}