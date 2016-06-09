# is.chebi
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Evaluates if a metabolite name is a ChEBI name
is.chebi<- function(metabolite){
  # Load ChEBI data
  chebi <- new.env()
  data("chebi", envir = chebi)
  # Return
  return(tolower(metabolite)%in%chebi$name)
}
