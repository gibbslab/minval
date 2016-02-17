# unbalanced
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Returns the unbalanced reactions from a set of stoichiometric reactions
unbalanced <- function(reaction, show.formulas = FALSE) {
  is.balanced <- function(reaction) {
    if (grepl("=", sub("(.*) <=> (.*)", "\\2", reaction))) {
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
    reactant <-
      unlist(mapply(
        function(coef, met) {
          rep(chebi.formula(.metname(met, rm.coef = TRUE)), coef)
        },
        coef = r_coef,
        met = reactant,
        SIMPLIFY = FALSE
      ))
    product <-
      unlist(mapply(
        function(coef, met) {
          rep(chebi.formula(.metname(met, rm.coef = TRUE)), coef)
        },
        coef = p_coef,
        met = product,
        SIMPLIFY = FALSE
      ))
    ! identical(.formula2matrix(reactant), .formula2matrix(product))
  }
  if (show.formulas == FALSE) {
    as.vector(sapply(reaction, is.balanced))
  } else {
    mb <- sapply(reaction, is.balanced)
    (cbind(reaction[mb], as.vector(sapply(reaction[mb], function(reaction) {
      toChEBI(reaction, formula = TRUE)
    }))))
  }
}
