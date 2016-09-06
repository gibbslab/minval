#' @importFrom stats complete.cases
#' @importFrom utils download.file read.delim2 read.table write.table
#' @importFrom XML saveXML
#' @importFrom XML xmlNode

# minval-internal
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Remove initial and final spaces for a metabolite
.remove.spaces <- function(metabolite){
  metabolite <- gsub("^[[:space:]]+","",metabolite)
  metabolite <- gsub("[[:space:]]+$","",metabolite)
  return(metabolite)
}

# Remove compartments for a metabolite
.remove.compartment <- function(met, rm.coef = FALSE) {
  met <- .remove.spaces(met)
  if (rm.coef == TRUE) {
    met <- .remove.coefficients(met)
  }
  gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$","",met)
}

# Remove coefficient for a metabolite
.remove.coefficients <- function(met){
  met <- gsub("^[[:digit:]]+[[:punct:]]*[[:digit:]]*[[:blank:]]+","",met)
  return(met)
}

.formula2matrix <- function(formula) {
  byatomtype <- unlist(regmatches(formula, gregexpr("([A-Z]{1}[a-z]?)([0-9]*)", formula)))
  atomtype <- sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', byatomtype)
  atomnumber <- as.numeric(regmatches(byatomtype, gregexpr('[0-9]+', byatomtype)))
  atomnumber[is.na(atomnumber)] <- 1
  tapply(atomnumber, atomtype, sum)
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

.get.right <- function(reaction){
  .remove.spaces(unlist(strsplit(unlist(strsplit(reaction,"[[:blank:]]*<?=>[[:blank:]]*"))[2],"[[:blank:]]+\\+[[:blank:]]+")))
}

.get.left <- function(reaction){
  .remove.spaces(unlist(strsplit(unlist(strsplit(reaction,"[[:blank:]]*<?=>[[:blank:]]*"))[1],"[[:blank:]]+\\+[[:blank:]]+")))
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

.write.tsv <- function(model,prefix){
 met <- matrix(as.vector(unlist(lapply(model$species, function(metabolite){unlist(metabolite)}))),ncol = 3,byrow = TRUE,dimnames = list(c(),c("abbreviation","name","compartment")))
 write.table(x = met,file = paste0(prefix,"_met.tsv"),row.names = FALSE)

 writeReaction <- function(reaction){
   compartment <- paste0(compartments(c(reaction[["reactants"]][["reactants"]],reaction[["products"]][["products"]])),collapse = ", ")
   reactants <- paste0(sapply(seq_along(reaction[["reactants"]][["reactants"]]), function(reactants){paste0("(",reaction[["reactants"]][["stoichiometry"]][reactants],") ",reaction[["reactants"]][["reactants"]][reactants])}),collapse = " + ")
   products <- paste0(sapply(seq_along(reaction[["products"]][["products"]]), function(products){paste0("(",reaction[["products"]][["stoichiometry"]][products],") ",reaction[["products"]][["products"]][products])}),collapse = " + ")
   id <- reaction[["id"]]
   reversible <- ifelse(reaction[["reversible"]],"reversible","irreversible")
   lb <- reaction[["parameters"]][["LOWER_BOUND"]]
   ub <- reaction[["parameters"]][["UPPER_BOUND"]]
   o <- reaction[["parameters"]][["OBJECTIVE_COEFFICIENT"]]
   rule <- reaction[["notes"]][["GPR"]]
   subsystem <- ""
   reaction <- paste(reactants,products,sep = ifelse(reaction[["reversible"]]," <==> "," --> "))
   reaction <- (c(id,id,reaction,reversible,compartment,lb,ub,o,rule,subsystem))
   return(reaction)
 }
 react <- matrix(unlist(sapply(model$reactions, writeReaction,simplify = FALSE)),ncol = 10,byrow = TRUE,dimnames = list(c(),c("abbreviation","name","equation","reversible","compartment","lowbnd","uppbnd","obj_coef","rule","subsystem")))
 write.table(x = react,file = paste0(prefix,"_react.tsv"),row.names = FALSE)
}

.sbmlCompatible <- function(metabolite,optimizedFor,type){
  compartment <- sapply(metabolite, compartments,USE.NAMES = FALSE)
  # Generando el ID
  if(optimizedFor == 'sybil' || optimizedFor == 'COBRA'){
    metabolite <- metabolites(metabolite,woCompartment = TRUE) 
  } else {
    metabolite <- metabolites(metabolite)
  }
  # Removiendo puntuaciones
  metabolite <- gsub("[[:blank:]]+","_",metabolite)
  metabolite <- paste0("M_",metabolite)
  metabolite <- gsub("[[:punct:]]+","_",metabolite)
  
  if (type == 's'){
    if(optimizedFor == 'sybil' || optimizedFor == 'COBRA'){
      metabolite <- paste0(metabolite,"[",compartment,"]")
      return(metabolite)
    } else {
      return (metabolite)
    }
  } else if (type == 'r'){
    if (optimizedFor =='RAVEN'){
      return(metabolite)
    } else{
      metabolite <- paste0(metabolite,"[",compartment,"]")
      return(metabolite)
    }
  }
}

.sbmlCompartment <- function(compartmentID, optimizedFor){
  if(optimizedFor == 'RAVEN'){
    return (paste0("C_",compartmentID))
  } else {
    return (compartmentID)
  }
}

.sbmlReaction <- function(reactionID, optimizedFor){
  if(optimizedFor == 'RAVEN'){
    return (paste0("R_",reactionID))
  } else {
    return (reactionID)
  }
}