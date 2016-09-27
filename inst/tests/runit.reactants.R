test.reactants<-function(){
  # Check output class
  checkTrue(is.vector(reactants("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
  checkTrue(is.list(reactants(glycolysis$REACTION)))
  # Check output values
  checkEquals(reactants("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ATP[m_a]","beta-D-Glucose[m_a]"))
}
