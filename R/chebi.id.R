# chebi.id
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Returns the ChEBI id asociated to a ChEBI name
chebi.id <- function(metabolite){
  # Load ChEBI data
  chebi <- new.env()
  data("chebi",package = "minval", envir = chebi)
  # Search in ChEBI names and returns the ChEBI id for a metabolite
  sapply(metabolite, function(metabolite){.safe.index(chebi$chebi[chebi$chebi$name%in%tolower(metabolite),1],1)})
}