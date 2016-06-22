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
  model <- .create.model()
  
  # Filling the model
  ## ID
  model$id <- sub("(.*)\\.(.*)$", "\\1", basename(infile))
  
  ## Compartments
  model$compartments <- lapply(compartments(data[,"EQUATION"]), .fill.compartment)
  
  ## Species
  metIDs <- metabolites(data[,"EQUATION"],woCompartment = FALSE,uniques = TRUE)
  model$species <- lapply(metIDs,.fill.species)
  
  ## Reactions
  model$reactions<-lapply(as.character(data[,"ID"]),.fill.reactions)
  
  # Writing model
  .write.xml(model,outfile)
}

