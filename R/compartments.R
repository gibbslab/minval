# compartments
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

compartments <- function(reactionList){
  # Extract metabolites
  metabolites <- metabolites(reactionList)
  # Extract compartments
  compartments <- unique(unlist(regmatches(metabolites, gregexpr("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$", metabolites))))
  # Remove brackets
  compartments <- gsub("\\[|\\]","",compartments)
  # Return compartments
  return(compartments)
}