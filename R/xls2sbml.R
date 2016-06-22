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
  for (compartment in compartments(data[,"EQUATION"])){
    model[["compartments"]][[length(model[["compartments"]])+1]] <- list(id=compartment,name=compartment)
  }
  
  ## Species
  for(metabolite in metabolites(data[,"EQUATION"],uniques = TRUE)){
    model[["species"]][[length(model[["species"]])+1]] <- list(id=met, name = metabolites(met,woCompartment = TRUE), compartment=compartments(met))
  }
  
  ## Reactions
  model$reactions<-lapply(as.character(data[,"ID"]),function(x){.fill.reactions(x,data)})
  
  # Writing model
  .write.xml(model,outfile)
}

