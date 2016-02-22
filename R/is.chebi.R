# is.chebi
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Evaluates if a metabolite name is a ChEBI name
is.chebi <- function(metabolite){
  # Load ChEBI data
  data("chebi", envir = environment())
  chebi <- chebi
  # Return
  if (is.na(table(chebi$name %in% tolower(metabolite))[2])){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
