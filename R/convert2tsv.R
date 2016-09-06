#' @export convert2tsv
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Write a TSV model for the 'sybil' R package.
#' @description This function converts a data.frame to a TSV model for the 'sybil' R package.
#' @details This function takes a data.frame as input and convert it to a valid sbmlR object, then the object is written into three \code{'.TSV'} output files (\code{'_react.tsv'}, \code{'_met.tsv'},\code{'_desc.tsv'}).
#' @param data A data.frame with the following mandatory colnames: \itemize{
#' \item \code{"ID":} A list of single character strings containing the reaction abbreviations, Entries in the field abbreviation are used as reaction ids, so they must be unique.
#' \item \code{"REACTION":} A set of stoichiometric reaction with the following format: \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.,
#' \item \code{"GPR":} A set of genes joined by boolean operators as AND or OR, rules may be nested by parenthesis. (optional: column can be empty),
#' \item \code{"LOWER.BOUND":} A list of numeric values containing the lower bounds of the reaction rates. 
#' If not set, zero is used for an irreversible reaction and 1000 for a reversible reaction. (optional: column can be empty),
#' \item \code{"UPPER.BOUND":} A list of numeric values containing the upper bounds of the reaction rates. 
#' If not set, 1000 is used by default. (optional: column can be empty),
#' \item \code{"OBJECTIVE":} A list of numeric values containing objective values for each reaction (optional: column can be empty). 
#' }
#' @param prefix A single character string in a writable path for three posible output '.TSV' files to be generated.
#' @return A set of three 'TSV' files in a valid format to the 'sybil' R package.
#' @examples  
#' \dontrun{
#' # Loading a CSV file
#' glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
#' 
#' # Data structure
#' head(glycolysis)
#' 
#' # Writing TSV files
#' convert2tsv(glycolysis,"glycolysis")
#' }
#' @keywords Convert TSV Metabolic Reconstruction
#' 
convert2tsv <- function(data, prefix){
  # Creating SBMLR model
  model <- convert2sbmlR(data,optimizedFor='sybil')
  # Function to write files
  write.tsv <- function(model,prefix){
    # Extracting metabolites
    met <- t(sapply(model$species,function(metabolite){c(metabolite[["name"]],metabolite[["compartment"]])}))
    met <- cbind(.sbmlCompatible(unique(met[,1])),unique(met[,1]),sapply(unique(met[,1]),function(metabolite){paste0(met[which(met[,1]==metabolite),2],collapse = ", ")},USE.NAMES = FALSE))
    colnames(met) <- c("abbreviation","name","compartment")
    # Writing metabolites file
    write.table(x = met,file = paste0(prefix,"_met.tsv"),row.names = FALSE,sep = "\t")
    # Extracting reaction data
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
    # Writing reaction file
    write.table(x = react,file = paste0(prefix,"_react.tsv"),row.names = FALSE,sep = "\t")
    # Extracting description information
    desc <- c(name = model$id,id = model$id)#, description = model$notes, compartment = paste0(sapply(model$compartments,function(compartment){compartment[["id"]]}),collapse = ", "), Nmetabolites = length(model$species), Nreactions=length(model$reactions))
    # Writing descriotion file
    write.table(t(desc),row.names = FALSE,file = paste0(prefix,"_desc.tsv"),sep = "\t")
  }
  # Main function
  write.tsv(model,prefix)
}
