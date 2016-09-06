#' @export mapReactions
#' @author Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Return reactions from a reference data.frame by a selected column
#' @description This function takes a set of values and return all associated data in a reference data.frame
#' 
#' @param reactionList A set of values to be mapped
#' @param referenceData A data.frame to be used as reference
#' @param by A data.frame column id or \code{'bool'}.
#' @param inverse A boolean value to define if inverted selection must be returned
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