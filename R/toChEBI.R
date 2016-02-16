# toChEBI
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Converts metabolite names to ChEBI ids in a stoichiometric reaction
toChEBI <- function(reaction, formula = FALSE) {
  reversible <- grepl("=", sub("(.*) <=> (.*)", "\\2", reaction))
  if (reversible) {
    reactant <-
      strsplit(sub("(.*) => (.*)", "\\1", reaction), " + ", fixed = TRUE)[[1]]
    product <-
      strsplit(sub("(.*) => (.*)", "\\2", reaction), " + ", fixed = TRUE)[[1]]
    
  } else {
    reactant <-
      strsplit(sub("(.*) <=> (.*)", "\\1", reaction), " + ", fixed = TRUE)[[1]]
    product <-
      strsplit(sub("(.*) <=> (.*)", "\\2", reaction), " + ", fixed = TRUE)[[1]]
  }
  r_coef <- as.numeric(sapply(reactant, .coeficients))
  r_coef[is.na(r_coef)] <- 1
  p_coef <- as.numeric(sapply(product, .coeficients))
  p_coef[is.na(p_coef)] <- 1
  
  if (formula == FALSE){
    reactant <- mapply(function(coef,met){paste(coef,chebi.id(.metname(met)),collapse =" ")}, coef=r_coef,met=reactant)
    product <- mapply(function(coef,met){paste(coef,chebi.id(.metname(met)),collapse =" ")}, coef=p_coef,met=product)
  } else {
    reactant <- mapply(function(coef,met){paste(coef,chebi.formula(.metname(met)),collapse =" ")}, coef=r_coef,met=reactant)
    product <- mapply(function(coef,met){paste(coef,chebi.formula(.metname(met)),collapse =" ")}, coef=p_coef,met=product)
  }
  
  if (reversible){
    paste(paste0(reactant,collapse = " + "),paste0(product, collapse = " + "), sep = " => ")
  } else {
    paste(paste0(reactant,collapse = " + "),paste0(product, collapse = " + "), sep = " <=> ")
  }
}

