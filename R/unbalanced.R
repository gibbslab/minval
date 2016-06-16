# unbalanced
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Returns the unbalanced reactions from a set of stoichiometric reactions
unbalanced <- function(reaction, show.formulas = FALSE) {
  ChEBIformulas <- toChEBI(reaction,formula = TRUE)
  r_met <- lapply(ChEBIformulas, .get.left)
  p_met <- lapply(ChEBIformulas, .get.right)
  r_met <- lapply(r_met, .atoms)
  p_met <- lapply(p_met, .atoms)
  r_met <- lapply(r_met, .formula2matrix)
  p_met <- lapply(p_met, .formula2matrix)
  balanced <- mapply(function(reactants,products){!identical(reactants,products)},reactants=r_met,products=p_met)
  if(length(grep("NA",ChEBIformulas))>0){
    warning("Some metabolites formulas were not found, mass unbalance was reported as NA",call. = FALSE)
  }
  balanced[grepl("NA",ChEBIformulas)] <- NA
  if (show.formulas==TRUE){
    ChEBIformulas[is.na(balanced)] <- "Some metabolites formulas were not found"
    balanced[grepl("NA",reaction)] <- FALSE
    balanced <- cbind(reaction=reaction[balanced==FALSE],formula=ChEBIformulas[balanced==FALSE])
  }
  return(balanced)
}
