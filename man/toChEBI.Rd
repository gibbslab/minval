\name{toChEBI}
\alias{toChEBI}
\title{ Translates compounds names to ChEBI ids in a stoichiometric reaction
}
\description{ This function translates compounds names to ChEBI ids or molecular formulas in a stoichiometric reaction}
\usage{
toChEBI(reaction, formula = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reaction}{A stoichiometric reaction}
  \item{formula}{ A logical value, if TRUE translates the compounds names to molecular formulas}
}
\examples{
## For a specific reaction
toChEBI("ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]")
# [1] "1 16761 + 1 18021 => 1 15422 + 1 15361"

toChEBI("ADP[c_n] + Phosphoenolpyruvate[c_n] => ATP[c_n] + Pyruvate[c_n]",formula = TRUE)
# [1] "1 C10H15N5O10P2 + 1 C3H4O6P => 1 C10H16N5O13P3 + 1 C3H3O3"

## For a set of reactions
# Loading a set of stoichiometric reactions from the glutamate/glutamine cycle
data("glugln")

# Translating to ChEBI ids
sapply(glugln[1:5],toChEBI,simplify = FALSE)
# $`alpha-D-glucose 1-phosphate[c_n] <=> alpha-D-Glucose 6-phosphate[c_n]`
# [1] "1 29042 <=> 1 17665"

# $`alpha-D-Glucose 6-phosphate[r_n] + water[r_n] => alpha-D-Glucose[r_n] + phosphate(3-)[r_n]`
# [1] "1 17665 + 1 15377 => 1 17925 + 1 18367"

# $`ATP[m_n] + alpha-D-Glucose[m_n] <=> ADP[m_n] + alpha-D-Glucose 6-phosphate[m_n]`
# [1] "1 15422 + 1 17925 <=> 1 16761 + 1 17665"

# $`ATP[m_n] + beta-D-Glucose[m_n] => ADP[m_n] + beta-D-Glucose 6-phosphate[m_n]`
# [1] "1 15422 + 1 15903 => 1 16761 + 1 17719"

# $`alpha-D-Glucose[m_n] + ADP[m_n] <=> alpha-D-Glucose 6-phosphate[m_n] + AMP[m_n]`
# [1] "1 17925 + 1 16761 <=> 1 17665 + 1 16027"

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
