#' @export validateSyntax
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Evaluate if a stoichiometric reaction has a valid syntax
#' @description  For a set of given stoichiometric reactions, this function makes the following syntactic evaluations for each reaction: \itemize{
#' \item Evaluates if the reaction contain more than one coefficient by metabolite
#' \item Evaluates if the reaction contain metabolite coefficients between parenthesis
#' \item Evaluates if the reaction contain arrow symbol between spaces
#' \item Evaluates if the reaction contain not allowed arrow symbols
#' \item Evaluates if the reaction contain metabolites name separated by a plus symbol between spaces
#' \item Evaluates if the reaction contain substituents separated of the metabolite names
#' }
#' @param reactionList A set of stoichiometric reaction with the following characteristics: \itemize{
#' \item Arrows symbols must be given in the form \code{'=>'} or \code{'<=>'}
#' \item Inverse arrow symbols \code{'<='} or other types as: \code{'-->'}, \code{'<==>'}, \code{'->'} will not be parsed and will lead to errors.
#' \item Arrow symbols and plus signs (\code{+}) must be surrounded by a space character
#' \item Stoichiometric coefficients must be surrounded by a space character and not by parentheses.
#' \item Each metabolite must have only one stoichiometric coefficient, substituents must be joined to metabolite name by a hyphen (\code{-}) symbol.
#' \item Exchange reactions have only one metabolite before arrow symbol
#' \item Compartments must be given between square brackets ([compartment]) joined at the end of metabolite name
#' }
#' Some examples of valid stoichiometric reactions are: \itemize{
#' \item \code{H2O[c] + Urea-1-Carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]}
#' \item \code{ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]}
#' \item \code{CO2[c] <=> }
#' }
#' @return  A boolean value \code{'TRUE'} if reaction has a valid syntax.
#' @examples 
#' # Evaluate the syntaxis for a single reaction
#' validateSyntax(reactionList = "ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]")
#' 
#' # Loading a set of stoichiometric reactions
#' glycolysis <- read.csv(system.file("extdata/glycolysisModel.csv",package = "minval"), sep='\t')
#' 
#' # Evaluating the syntaxis for a set of stoichiometric reactions
#' validateSyntax(reactionList = glycolysis$REACTION)

validateSyntax <- function(reactionList){
  # Convert to a vector
  reactionList <- as.vector(reactionList)
  # Empty vector
  valid.syntax <- NULL
  # Coefficient validation
  valid.syntax <- c(valid.syntax,grepl(pattern = "([[:digit:]][[:space:]][[:digit:]][[:space:]])+",x = reactionList))
  valid.syntax <- c(valid.syntax,grepl(pattern = "(\\([[:digit:]]+\\)[[:space:]]+)",x = reactionList))
  # Directionality validation
  valid.syntax <- c(valid.syntax, (!grepl(pattern = "<?=>?",x = reactionList)))
  valid.syntax <- c(valid.syntax, (!grepl(pattern = "[[:graph:]]+[[:space:]]+<?=>[[:print:]]*",x = reactionList)))
  valid.syntax <- c(valid.syntax,grepl(pattern = "(<)?\\-(\\-)?>",x = reactionList))
  valid.syntax <- c(valid.syntax,(grepl(pattern = "<-(>)?",x = reactionList) | grepl(pattern = "<=[[:space:]]+",x = reactionList)))
  # Metabolite names validation
  valid.syntax <- c(valid.syntax,grepl(pattern = "[[:alnum:]]+\\+[[:alnum:]]+",x = reactionList))
  # Blank spaces validation
  valid.syntax <- c(valid.syntax, (!grepl(pattern = "[[:space:]]",x = reactionList)))
  # Validating metabolite name
  valid.syntax <- c(valid.syntax, grepl(pattern = "[[:space:]]\\-[[:alnum:]]",x = reactionList))
  # Metabolite composition
  nMetabolites <- lengths(lapply(reactionList,function(reaction){metabolites(reaction)}))
  valid.syntax <- c(valid.syntax, (nMetabolites < 1))
  valid.syntax <- c(valid.syntax, unlist(lapply(reactionList, function(reaction){
    if(length(metabolites(reaction))>0){
      left <- metabolites(unlist(getLeft(reaction)))
      left <- length(left[!(is.null(left) || is.na(left))])
      right <- metabolites(unlist(getRight(reaction)))
      right <- length(right[!(is.null(right) || is.na(right))])
      return(!(left >= 1 & right > 0 | left == 1 & right == 0))
    } else {
      return(FALSE)
    }
  })))
  # Warnings!
  valid.syntax <- matrix(valid.syntax,ncol = 11)
  sapply(which(valid.syntax[,1]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid coefficients. Metabolites should have just one coefficient."),call. = FALSE)})
  sapply(which(valid.syntax[,2]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid coefficients. Coefficients should have not parentheses."),call. = FALSE)})
  sapply(which(valid.syntax[,3]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Not arrow detected."),call. = FALSE)})
  sapply(which(valid.syntax[,4]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Arrow symbols should be separated of metabolites by a blank space."),call. = FALSE)})
  sapply(which(valid.syntax[,5]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Please use <=> or => instead of <-> or -> or -->."),call. = FALSE)})
  sapply(which(valid.syntax[,6]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid directionality symbols. Reverse symbols <= or <- are not allowed"),call. = FALSE)})
  sapply(which(valid.syntax[,7]==TRUE),function(x){warning(paste0("Reaction ",x,": Metabolites names should be separated by a plus symbol between spaces."),call. = FALSE)})
  sapply(which(valid.syntax[,8]==TRUE),function(x){warning(paste0("Reaction ",x,": No blank spaces detected between metabolites."),call. = FALSE)})
  sapply(which(valid.syntax[,9]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid metabolite name. Substituent should not be separated of the metabolite name."),call. = FALSE)})
  sapply(which(valid.syntax[,10]==TRUE),function(x){warning(paste0("Reaction ",x,": Invalid metabolite composition. There are not metabolites."),call. = FALSE)})
  sapply(which(valid.syntax[,11]==TRUE),function(x){warning(paste0("Reaction ",x,": Not valid syntax. Exchange reactions should have only one metabolite before arrow symbol"),call. = FALSE)})
  # Return
  if(any(valid.syntax == TRUE)){
    stop("Please check that:\n* Arrows and plus signs are surrounded by a 'space character'\n* Stoichiometric coefficients are surrounded by spaces\n* Arrows to be in the form '=>' or '<=>'\n* Exchange reactions have only one metabolite before arrow symbol",call. = FALSE)
  } else {
    return(sapply(seq_along(reactionList), function(reaction){all(valid.syntax[reaction,]==FALSE)}))
  }
}
