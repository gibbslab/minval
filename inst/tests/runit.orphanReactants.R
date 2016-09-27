test.orphanReactants<-function(){
  # Check output class
  checkTrue(is.vector(orphanReactants("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(orphanReactants("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ATP[m_a]","beta-D-Glucose[m_a]"))
}
