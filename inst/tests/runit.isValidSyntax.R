test.isValidSyntax<-function(){
  # Check output class
  checkTrue(is.vector(isValidSyntax("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
  checkTrue(is.vector(isValidSyntax(glycolysis$REACTION)))
  
  # Check output values
  checkEquals(isValidSyntax("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),TRUE)
  

}