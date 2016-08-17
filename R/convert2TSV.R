# convert2TSV
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

convert2TSV <- function(data, prefix){
  # Creating SBMLR model
  model <- convert2SBMLR(data)
  # Function to write files
  write.tsv <- function(model,prefix){
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
  write.tsv(model,prefix)
}
