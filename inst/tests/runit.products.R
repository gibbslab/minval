test.products<-function(){
  # Check output class
  checkTrue(is.vector(products("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
  checkTrue(is.list(products(glycolysis$REACTION)))
  
  # Check output values
  checkEquals(products("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ADP[m_a]","beta-D-Glucose 6-phosphate[m_a]"))
}
