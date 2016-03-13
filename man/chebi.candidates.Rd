\name{chebi.candidates}
\alias{chebi.candidates}
\title{Returns the possible ChEBI names based on compound synonyms
}
\description{
This function takes a compound or a set of compounds names as an input and returns a list of possible CheBI names.
}
\usage{
chebi.candidates(metabolite)
}
\arguments{
  \item{metabolite}{A compound or a set of compounds names}
}
\examples{
## For a specific metabolite
chebi.candidates("veg")
# [1] "lutein"
# [2] "EC 2.7.10.1 (receptor protein-tyrosine kinase) inhibitor"
# [3] "vascular endothelial growth factor receptor antagonist"
# [4] "paliperidone"
# [5] "paliperidone palmitate"

## For a set of metabolites
# Loading a set of stoichiometric reactions from the glutamate/glutamine cycle
data("glugln")

# Extracting metabolites names
mets<-metabolites(glugln[1:2],woCompartment = TRUE)

# Searching for candidates
sapply(mets,chebi.candidates)

# $`alpha-D-glucose 1-phosphate`
# [1] "alpha-D-glucose 1-phosphate"
# [2] "alpha-D-glucose 1-phosphate(2-)"

# $`alpha-D-Glucose 6-phosphate`
# [1] "alpha-D-glucose 6-phosphate"
# [2] "D-glucopyranose 6-phosphate"
# [3] "alpha-D-glucose 6-phosphate(2-)"

# $water
#  [1] "water"
#  [2] "phenol"
#  [3] "ammonia"
#  [4] "hydrogen sulfide"
#  [5] ...
#[165] "spectinomycin hydrochloride hydrate"

#$`alpha-D-Glucose`
# [1] "dTDP-3-amino-3,4,6-trideoxy-alpha-D-glucose"
# [2] "dTDP-alpha-D-glucose"
# [3] "D-glucopyranose 1-phosphate"
# [4] "dTDP-4-dehydro-6-deoxy-alpha-D-glucose"
# [5] "alpha-D-glucose 6-phosphate"
#[79] "lipid A-F"

#$`phosphate(3-)`
#character(0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
