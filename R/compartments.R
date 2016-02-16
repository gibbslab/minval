# compartments
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

compartments <- function(metabolites){
  # Compartment main function
  compartment<- function(met){
    # Extract the letter between the [ ]
    strsplit(strsplit(met,fixed=TRUE,split="[")[[1]][length(strsplit(met,fixed=TRUE,split="[")[[1]])],fixed = TRUE,split = "]")[[1]]
  }
  # Extract the compartments for a set of stoichiometric reaction
  unique(sapply(metabolites,compartment))
}
