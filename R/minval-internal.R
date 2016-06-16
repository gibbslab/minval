# minval-internal
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

.remove.compartment <- function(met, rm.coef = FALSE) {
  met <- gsub("^[[:blank:]]*","",met)
  met <- gsub("[[:blank:]]*$","",met)
  if (rm.coef == TRUE) {
    met <- gsub("^[[:digit:]][[:graph:]]*[[:blank:]]","",met)
  }
    gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$","",met)
}

.formula2matrix <- function(formula) {
  byatomtype <- unlist(regmatches(formula, gregexpr("([[:alpha:]]{1}[[:alpha:]]*)([0-9]*)", formula)))
  atomtype <- sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', byatomtype)
  atomnumber <- as.numeric(regmatches(byatomtype, gregexpr('[0-9]+', byatomtype)))
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
  met <- as.numeric(met)
  return(met)
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

.get.right <- function(reaction){
  unlist(strsplit(unlist(strsplit(reaction,"[[:blank:]]+<?=>[[:blank:]]+"))[2],"[[:blank:]]+\\+[[:blank:]]+"))
}

.get.left <- function(reaction){
  unlist(strsplit(unlist(strsplit(reaction,"[[:blank:]]+<?=>[[:blank:]]+"))[1],"[[:blank:]]+\\+[[:blank:]]+"))
}

.join.cm <- function(coefficient,metabolite){
  joined <- paste(mapply(function(coefficient,metabolite){paste(coefficient,metabolite,collapse =" ")},coefficient=coefficient,metabolite=metabolite),collapse = " + ")
  return(joined)
}

.join.reaction <- function(reactant,reversibility,product){
  if (reversibility) {
    reaction <- paste(reactant,product,sep = " <=> ")
  } else {
    reaction <- paste(reactant,product,sep = " => ")
  }
}