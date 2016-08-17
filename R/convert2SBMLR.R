# convert2SBMLR
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

convert2SBMLR <- function(data){
  # Import data.frame
  data <- as.data.frame.array(data)
  # Validate colnames
  .validate.xls(data)
  # Remove comments
  data <- .remove.comments(data)
  # Verify syntaxis
  data <- data[is.validSyntax(as.vector(data[,"REACTION"])),]
  # Remove arrows
  data[,"REACTION"] <- gsub("<?->","-",data[,"REACTION"])
  # Creating the model
  model <- .create.model()
  # Filling the model
  ## ID
  model$id <- "model"
  ## Compartments
  model$compartments <- lapply(compartments(data[,"REACTION"]),function(compartment){list(id=compartment,name=compartment)}) 
  ## Species
  model$species <- lapply(metabolites(data[,"REACTION"],uniques = TRUE),function(met){list(id=met, name = metabolites(met,woCompartment = TRUE), compartment=compartments(met))})
  ## Reactions
  model$reactions<-lapply(as.character(data[,"ID"]),function(x){.fill.reactions(x,data)})
  # Return
  return(model)
}