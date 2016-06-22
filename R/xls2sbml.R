# xls2sbml
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

xls2sbml<-function(infile,outfile){
  
  # Getting data
  data <- as.data.frame.array(read.xls(infile , sheet = 1))
  .validate.xls(data)
  data <- .remove.comments(data)
  data[,"EQUATION"] <- gsub("<?->","-",data[,"EQUATION"])
  
  # Creating the model
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
  
  # Filling the model
  ## ID
  model$id <- sub("(.*)\\.(.*)$", "\\1", basename(infile))
  
  ## Compartments
  .fill.compartment <- function(compartment){
    model[[3]][[length(model[[3]])+1]] <- list(id=compartment,name=compartment)
  }
  model$compartments <- lapply(compartments(data[,"EQUATION"]), .fill.compartment)
  
  ## Species
  .fill.species<- function(met){ 
    model[[4]][[length(model[[4]])+1]] <- list(id=met, name = metabolites(met,woCompartment = TRUE), compartment=compartments(met)
    )
  }
  model$species <- lapply(metabolites(data[,"EQUATION"],woCompartment = FALSE,uniques = TRUE),.fill.species)
  
  ## Reactions
  model$reactions<-lapply(as.character(data[,"ID"]),function(x){.fill.reactions(x,data)})
  
  # Writing model
  .write.xml(model,outfile)
}

