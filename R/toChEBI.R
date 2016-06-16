# toChEBI
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Converts metabolite names to ChEBI ids in a stoichiometric reaction
toChEBI <- function(reaction, formula = FALSE) {
  chebi <- new.env()
  data("chebi",package = "minval", envir = chebi)
  .chebi <- function(metabolite,get){
    data <- structure(chebi$chebi[,get],names=chebi$chebi$name)
    return (as.vector(data[tolower(metabolite)]))
  }
  # Evaluates reversibility
  reversible <- grepl("<=>",reaction)
  # Extract metabolites
  r_met <- lapply(reaction, .get.left)
  p_met <- lapply(reaction, .get.right)
  # Extract coefficient
  r_coef <- lapply(r_met, .coefficients)
  p_coef <- lapply(p_met, .coefficients)
  # Remove coefficient and compartments
  r_met <- lapply(r_met, function(reaction){metabolites(reaction,woCompartment = TRUE,uniques = FALSE)})
  p_met <- lapply(p_met, function(reaction){metabolites(reaction,woCompartment = TRUE,uniques = FALSE)})
  # Find associated data
  if (formula == FALSE){
    r_met <- lapply(r_met, function(metabolites){as.vector(.chebi(metabolites,get="id"))})
    p_met <- lapply(p_met, function(metabolites){as.vector(.chebi(metabolites,get="id"))})
  } else {
    r_met <- lapply(r_met, function(metabolites){as.vector(.chebi(metabolites,get="formula"))})
    p_met <- lapply(p_met, function(metabolites){as.vector(.chebi(metabolites,get="formula"))})
  }
  # Join metabolites and coefficient
  r_met <- mapply(.join.cm,coefficient=r_coef,metabolite=r_met,USE.NAMES = FALSE)
  p_met <- mapply(.join.cm,coefficient=p_coef,metabolite=p_met,USE.NAMES = FALSE)
  # Join converted reaction
  reaction <- mapply(.join.reaction,reactant=r_met,reversibility=reversible,product=p_met,USE.NAMES = FALSE)
  # Return converted reaction
  return(reaction)
}

