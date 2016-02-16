# is.chebi
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Evaluates if a metabolite name is a ChEBI name
is.chebi <- function(metabolite){
  # Load ChEBI data
  data("chebi", envir = environment())
  chebi <- chebi
  chebi.match<-chebi[match(tolower(metabolite),tolower(chebi$name)),]
  if (is.na(chebi.match[1])){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
