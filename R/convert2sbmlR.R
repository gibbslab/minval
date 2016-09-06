#' @export convert2sbmlR
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Convert a data.frame data to a SBMLR object
#' @description This function converts a data.frame to a SBML-like R list of lists core object of class SBMLR.
#' The Systems Biology Markup Language (SBML) is a representation format, based on XML, for communicating and storing computational models of biological processes.
#' More Info: Encyclopedia of Systems Biology Dubitzky, W., Wolkenhauer, O., Yokota, H., Cho, K.-H. (Eds.) SBML, pp2057-2062 Springer 2013.
#' @param data A data.frame with the following mandatory colnames: \itemize{
#' \item \code{"ID":} A list of single character strings containing the reaction abbreviations, Entries in the field abbreviation are used as reaction ids, so they must be unique.
#' \item \code{"REACTION":} A set of stoichiometric reaction with the following format: \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}".
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.,
#' \item \code{"GPR":} A set of genes joined by boolean operators as AND or OR, rules may be nested by parenthesis. (optional: column can be empty),
#' \item \code{"LOWER.BOUND":} A list of numeric values containing the lower bounds of the reaction rates.
#' If not set, zero is used for an irreversible reaction and 1000 for a reversible reaction. (optional: column can be empty),
#' \item \code{"UPPER.BOUND":} A list of numeric values containing the upper bounds of the reaction rates.
#' If not set, 1000 is used by default. (optional: column can be empty),
#' \item \code{"OBJECTIVE":} A list of numeric values containing objective values for each reaction (optional: column can be empty).
#' }
#' @param optimizedFor A character string specifying the toolbox for which the SBML file must be optimized; must be one of \code{'sybil'}, \code{'RAVEN'} or \code{'COBRA'}.
#' @return A SBML-like R list of lists core object of class SBMLR
convert2sbmlR <- function(data,optimizedFor){
  # Import data.frame
  data <- as.data.frame.array(data)
  # Validate colnames
  validateData <- function(data){
    names <- colnames(data)
    if (length(grep("^ID$",names,ignore.case = TRUE))==0){stop("Reaction ID's not found")}
    if (!identical(data[,"ID"],unique(data[,"ID"]))){stop("Reaction ID's must be unique")}
    if (length(grep("^REACTION$",names,ignore.case = TRUE))==0){stop("REACTIONS not found")}
    if (length(grep("^GPR$",names,ignore.case = TRUE))==0){stop("GPR not found")}
    if (length(grep("^LOWER.BOUND$",names,ignore.case = TRUE))==0){stop("LB not found")}
    if (length(grep("^UPPER.BOUND$",names,ignore.case = TRUE))==0){stop("UB not found")}
    if (length(grep("^OBJECTIVE$",names,ignore.case = TRUE))==0){stop("OBJECTIVE not found")}
  }
  validateData(data)
  # Remove comments
  removeComments <- function(data){
    if (length(grep("#",data[,1]))>0){
      data <- data[!data[,1]=="#"]
    }
    return (data)
  }
  data <- removeComments(data)
  # Verify syntaxis
  validSyntax <- isValidSyntax(data[,"REACTION"])
  data <- data[validSyntax,]
  # Remove arrows
  data[,"REACTION"] <- gsub("<?->","-",data[,"REACTION"])
  # Creating the model
  createModel <- function(){
    model <- list(
      id = "",
      notes = c(""),
      compartments = list(),
      species = list(),
      reactions = list(),
      globalParameters = list(),
      rules = list()
    )
    model <-structure(model,class ="SBMLR")
    return(model)
  }
  model <- createModel()
  # Filling the model
  ## ID
  model$id <- "model"
  ## Compartments
  model$compartments <- lapply(compartments(data[,"REACTION"]),function(compartment){list(id=.sbmlCompartment(compartment,optimizedFor),name=compartment)})
  ## Species
  model$species <- lapply(metabolites(data[,"REACTION"],uniques = TRUE),function(met){list(id=.sbmlCompatible(met,optimizedFor,'s'), name = metabolites(met,woCompartment = TRUE), compartment=.sbmlCompartment(compartments(met),optimizedFor))})
  ## Reactions
  fillReactions <- function(rxnid,data){
    LB = ifelse(is.na(data[data[,"ID"]%in%rxnid,"LOWER.BOUND"]),ifelse(grepl("<=>",data[data[,"ID"]%in%rxnid,"REACTION"]),-1000,0),data[data[,"ID"]%in%rxnid,"LOWER.BOUND"])
    UB = ifelse(is.na(data[data[,"ID"]%in%rxnid,"UPPER.BOUND"]),1000,data[data[,"ID"]%in%rxnid,"UPPER.BOUND"])
    rev <- grepl("<=>",data[data[,"ID"]%in%rxnid,"REACTION"])
    left <- .get.left(data[data[,"ID"]%in%rxnid,"REACTION"])
    right <- .get.right(data[data[,"ID"]%in%rxnid,"REACTION"])
    gpr <- data[data[,"ID"]%in%rxnid,"GPR"]
    genes <- unlist(strsplit(gsub("[(and|or)]","",gpr),"[[:blank:]]+"))
    reac<-list(id=.sbmlReaction(as.vector(rxnid),optimizedFor),
               reversible=rev,
               reactants=list(reactants=.sbmlCompatible(left,optimizedFor,'r'),
                              stoichiometry=.coefficients(left)),
               products=list(products=ifelse(is.na(right),gsub("\\[[[:graph:]]+\\]","\\[b\\]",.sbmlCompatible(left,optimizedFor,'r')),.sbmlCompatible(right,optimizedFor,"r")),stoichiometry=.coefficients(right)),
               parameters=c(LOWER_BOUND = LB,
                            UPPER_BOUND = UB,
                            OBJECTIVE_COEFFICIENT = ifelse(is.na(data[data[,"ID"]%in%rxnid,"OBJECTIVE"]),0,data[data[,"ID"]%in%rxnid,"OBJECTIVE"]),
                            FLUX_VALUE = 0),
               mathmlLaw = XML::xmlNode("ci","FLUX_VALUE"),
               strlaw = "FLUX_VALUE",
               notes=list(GPR=gpr,
                          GENE=genes))
    return(reac)
  }
  model$reactions<-lapply(as.character(data[,"ID"]),function(x){fillReactions(x,data)})
  # Return
  return(model)
}
