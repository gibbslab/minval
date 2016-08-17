# convert2SBMLR
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

convert2SBMLR <- function(data){
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
  validSyntax <- is.validSyntax(data[,"REACTION"])
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
  model$compartments <- lapply(compartments(data[,"REACTION"]),function(compartment){list(id=compartment,name=compartment)}) 
  ## Species
  model$species <- lapply(metabolites(data[,"REACTION"],uniques = TRUE),function(met){list(id=met, name = metabolites(met,woCompartment = TRUE), compartment=compartments(met))})
  ## Reactions
  fillReactions <- function(rxnid,data){
    LB = ifelse(is.na(data[data[,"ID"]%in%rxnid,"LOWER.BOUND"]),ifelse(grepl("<=>",data[data[,"ID"]%in%rxnid,"REACTION"]),-1000,0),data[data[,"ID"]%in%rxnid,"LOWER.BOUND"])
    UB = ifelse(is.na(data[data[,"ID"]%in%rxnid,"UPPER.BOUND"]),1000,data[data[,"ID"]%in%rxnid,"UPPER.BOUND"])
    rev <- grepl("<=>",data[data[,"ID"]%in%rxnid,"REACTION"])
    left <- .get.left(data[data[,"ID"]%in%rxnid,"REACTION"])
    right <- .get.right(data[data[,"ID"]%in%rxnid,"REACTION"])
    gpr <- data[data[,"ID"]%in%rxnid,"GPR"]
    genes <- unlist(strsplit(gsub("[(and|or)]","",gpr),"[[:blank:]]+"))
    reac<-list(id=as.vector(rxnid), 
               reversible=rev,
               reactants=list(reactants=metabolites(left),stoichiometry=.coefficients(left)),
               products=list(products=ifelse(is.na(right),paste0(metabolites(left,woCompartment=TRUE),"[b]"),metabolites(right)),stoichiometry=.coefficients(right)),
               parameters=c(LOWER_BOUND = LB,
                            UPPER_BOUND = UB,
                            OBJECTIVE_COEFFICIENT = ifelse(is.na(data[data[,"ID"]%in%rxnid,"OBJECTIVE"]),0,data[data[,"ID"]%in%rxnid,"OBJECTIVE"]),
                            FLUX_VALUE = 0),
               mathmlLaw = xmlNode("ci","FLUX_VALUE"),
               strlaw = "FLUX_VALUE",
               notes=list(GPR=gpr,
                          GENE=genes))
    return(reac)
  }
  model$reactions<-lapply(as.character(data[,"ID"]),function(x){fillReactions(x,data)})
  # Return
  return(model)
}