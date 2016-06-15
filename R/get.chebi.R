# get.chebi
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

get.chebi <- function(metabolite,get="id"){
  chebi <- new.env()
  data("chebi",package = "minval", envir = chebi)
  get <- match.arg(get,c("id","formula","candidates"))
  if(get=="candidates"){
    candidates <- sapply(metabolite, function(metabolite){chebi$chebi$name[grep(metabolite,chebi$chebi$candidates,ignore.case = TRUE)]},simplify = FALSE)
    return(candidates)
  } else {
    name <- tolower(metabolite)
    data <- structure(chebi$chebi[,get],names=chebi$chebi$name)
    data <- structure(data[metabolite],names=metabolite)
    return (data)
  }
}
