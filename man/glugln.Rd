\name{glugln}
\alias{glugln}
\docType{data}
\title{Stoichiometric reactions from the reconstructon of the glutamate/glutamine cycle
%%   ~~ data name/kind ... ~~
}
\description{
This dataset contains the stoichiometric reactions from the reconstructon of the glutamate/glutamine cycle. It was presented in a workshop of the "Latin-American School on glial cells in the diseased brain" to understand how glutamate can be a neurotoxic agent through glutamate-induced excitotoxicity in neurons and astrocytes. All the reactions have the following format: L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
}
\usage{data("glugln")}
\format{
 An array with 128 compartmentalized stoichiometric reactions from the metabolic reconstruction of Glutamate-Glutamine cycle.
}
\source{
Vega-Vela N, Jimenez CM, Barreto G and Gonzalez J (2015). Metabolic Reconstruction of Glutamate-Glutamine Cycling: A Flux Balance Approach. Front. Cell. Neurosci. Conference Abstract: Latin-American School on glial cells in the diseased brain (IBRO). doi: 10.3389/conf.fncel.2015.35.00008
}
\examples{
data("glugln")
head(glugln)
# [1] "D-Glucose 1-phosphate[c_n] <=> alpha-D-Glucose 6-phosphate[c_n]"
# [2] "alpha-D-Glucose 6-phosphate[r_n] + H2O[r_n] => alpha-D-Glucose[r_n] + Orthophosphate[r_n]"
# [3] "ATP[m_n] + alpha-D-Glucose[m_n] <=> ADP[m_n] + alpha-D-Glucose 6-phosphate[m_n]"
# [4] "ATP[m_n] + beta-D-Glucose[m_n] => ADP[m_n] + beta-D-Glucose 6-phosphate[m_n]"
# [5] "alpha-D-Glucose[m_n] + ADP[m_n] <=> alpha-D-Glucose 6-phosphate[m_n] + AMP[m_n]"
# [6] "beta-D-Glucose[m_n] + ADP[m_n] <=> beta-D-Glucose 6-phosphate[m_n] + AMP[m_n]"
}
\keyword{datasets}
