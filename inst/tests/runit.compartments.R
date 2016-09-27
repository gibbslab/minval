test.compartments<- function(){
  # Check output class
  checkTrue(is.vector(compartments("A[c]")))
  # Check values
  checkEquals(compartments(c("A[c]","B[protein]A[m_pro]","C[r]")),c("c","m_pro","r"))
}