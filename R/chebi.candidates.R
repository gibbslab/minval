# chebi.candidates
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

## Returns the possible ChEBI names based on metabolite synonyms
chebi.candidates <- function(metabolite) {
  # Load ChEBI data
  chebi <- new.env()
  data("chebi",package = "minval", envir = chebi)
  # Search in metabolite synonyms the metabolite name
  sapply(metabolite, function(metabolite){chebi$chebi$name[grep(metabolite, chebi$chebi$synonyms, ignore.case = TRUE)]})
}
