# xls2tsv
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

xls2tsv <- function(infile, prefix){
  # Reading data
  data <- gdata::read.xls(infile,sheet = 1)
  # Creating SBMLR model
  model <- convert2sbmlR(data)
  # Function to write files
  write.tsv <- function(model,prefix){
    met <- t(sapply(model$species,function(metabolite){c(metabolite[["name"]],metabolite[["compartment"]])}))
    met <- cbind(unique(met[,1]),unique(met[,1]),sapply(unique(met[,1]),function(metabolite){paste0(met[which(met[,1]==metabolite),2],collapse = ", ")},USE.NAMES = FALSE))
    colnames(met) <- c("abbreviation","name","compartment")
    write.table(x = met,file = paste0(prefix,"_met.tsv"),row.names = FALSE,sep = "\t")
    
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
    write.table(x = react,file = paste0(prefix,"_react.tsv"),row.names = FALSE,sep = "\t")
    desc <- c(name = model$id,id = model$id)#, description = model$notes, compartment = paste0(sapply(model$compartments,function(compartment){compartment[["id"]]}),collapse = ", "), Nmetabolites = length(model$species), Nreactions=length(model$reactions))
    write.table(t(desc),row.names = FALSE,file = paste0(prefix,"_desc.tsv"),sep = "\t")
  }
  write.tsv(model,prefix)
}
