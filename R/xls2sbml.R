# xls2sbml
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

xls2sbml<-function(infile,outfile){
  
  # Getting data
  data <- as.data.frame.array(read.xls(infile , sheet = 1))
  .validate.xls(data)
  data <- .remove.comments(data)
  data[,"REACTION"] <- gsub("<?->","-",data[,"REACTION"])
  
  # Creating the model
  model <- .create.model()
  
  # Filling the model
  ## ID
  model$id <- sub("(.*)\\.(.*)$", "\\1", basename(infile))
  
  ## Compartments
  model$compartments <- lapply(compartments(data[,"REACTION"]),function(compartment){list(id=compartment,name=compartment)}) 
  
  ## Species
  model$species <- lapply(metabolites(data[,"REACTION"],uniques = TRUE),function(met){list(id=met, name = metabolites(met,woCompartment = TRUE), compartment=compartments(met))})
  
  ## Reactions
  model$reactions<-lapply(as.character(data[,"ID"]),function(x){.fill.reactions(x,data)})
  
  # Writing model
  .write.xml(model,outfile)
}

