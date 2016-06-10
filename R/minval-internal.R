# minval-internal
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

.metname <- function(met, rm.coef = FALSE) {
  met <- gsub("^[[:blank:]]*","",met)
  met <- gsub("[[:blank:]]*$","",met)
  if (rm.coef == TRUE) {
    met <- gsub("^[[:digit:]][[:graph:]]*[[:blank:]]","",met)
  }
    gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$","",met)
}

.formula2matrix <- function(formula) {
  byatomtype <-
    unlist(regmatches(formula, gregexpr("([A-Z]{1}[a-z]?)([0-9]*)", formula)))
  atomtype <- sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', byatomtype)
  atomnumber <-
    as.numeric(regmatches(byatomtype, gregexpr('[0-9]+', byatomtype)))
  atomnumber[is.na(atomnumber)] <- 1
  tapply(atomnumber, atomtype, sum)
}

.remove.coefficients <- function(met){
  met <- gsub("^[[:digit:]]+[[:punct:]]*[[:digit:]]*[[:blank:]]+","",met)
  return(met)
}

.coefficients <- function(met) {
  met <- regmatches(met, gregexpr('^[[:digit:]][[:punct:]]*[[:digit:]]*[[:blank:]]+', met))
  met[lengths(met)==0] <- 1
  met <- gsub("[[:blank:]]*$","",met)
  return(met)
}

.safe.index <- function(df, n){
  tryCatch(df[n], error = function(e)return(rep(NA, nrow(df))))
}

.atoms <- function(metabolites) {
  coef <- as.numeric(sapply(metabolites, .coefficients))
  formula <- metabolites(metabolites)
  unlist(mapply(function(coef, formula){rep(formula,coef)}, coef = coef, formula = formula,SIMPLIFY = FALSE))
}

.remove.spaces <- function(metabolite){
  metabolite <- gsub("^[[:space:]]","",metabolite)
  metabolite <- gsub("[[:space:]]$","",metabolite)
  return(metabolite)
}
