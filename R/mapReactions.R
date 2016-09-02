#' @export mapReactions
#' @author Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Return reactions of a reference data by a selected column
#' @description
#' 
#' @param reactionList 
#' @param referenceData 
#' @param by 
#' @param inverse 
#' 
#' @examples 
#' # Loading a CSV file
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#'
#' # Data structure
#' head(glycolysis)
#' 
#' # Mapping reactions
#' mapReactions(reactionList = isValidSyntax(glycolysis$REACTION), referenceData = glycolysis, by = "bool")
#' 
#' # Mapping inverse
#' mapReactions(reactionList = isValidSyntax(glycolysis$REACTION), referenceData = glycolysis, by = "bool", inverse = TRUE)

mapReactions <- function(reactionList, referenceData, by, inverse=FALSE ){
  if(!is.null(dim(referenceData)) && is.null(dim(reactionList))){
    referenceData <- as.data.frame.array(referenceData)
    if(by=="bool"){
      if (inverse == FALSE){
        return(referenceData[reactionList,])
      } else {
        return(referenceData[!reactionList,])
      }
    } else if(by%in%colnames(referenceData)){
      if(inverse == FALSE){
        return(referenceData[referenceData[,match(by,colnames(referenceData))]%in%reactionList,])
      } else {
        return(referenceData[!referenceData[,match(by,colnames(referenceData))]%in%reactionList,])
      }
    }
  }
}