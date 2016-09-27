test.stoichiometricMatrix<- function(){
  # Check output class
  checkTrue(is.matrix(stoichiometricMatrix("A[c] + B[c] => D[c]")))
  # Check values
  checkEquals(dim(stoichiometricMatrix("A[c] + B[c] => D[c]")),c(3,1))
}