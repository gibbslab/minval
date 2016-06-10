# chebi.formula
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

## Returns the molecular formula associated to a ChEBI metabolite name
chebi.formula <- function(metabolite){
  # Load ChEBI data
  chebi <- new.env()
  data("chebi",package = "minval", envir = chebi)
  # Search in ChEBI database for a molecular formula based in metabolite name
  sapply(metabolite, function(metabolite){.safe.index(chebi$chebi[chebi$chebi$name%in%tolower(metabolite),4],1)})
}
