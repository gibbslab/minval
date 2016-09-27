test.isBalanced<-function(){
  chemicalData <- read.csv2(system.file("extdata", "chemData.csv", package = "minval"))
  glycolysis <- read.csv2(system.file("extdata", "glycolysisKEGG.csv", package = "minval"))
  
  
  # Check output class
  checkTrue(is.vector(isBalanced(reactionList = glycolysis$REACTION,referenceData = chemicalData,mFormula = "FORMULA",ids = "NAME")))
  
  # Check output values
  checkEquals(isBalanced(reactionList = "alpha-D-Glucose[c] + Orthophosphate[c] => alpha-D-Glucose 6-phosphate[c] + H2O[c]",referenceData = chemicalData,mFormula = "FORMULA",ids = "NAME"),TRUE)
  checkEquals(isBalanced(reactionList = "ATP[c] + Phosphoenolpyruvate[c] => Pyruvate[c] + ADP[c]",referenceData = chemicalData,mFormula = "FORMULA",ids = "NAME"),FALSE)
  checkEquals(isBalanced(reactionList = glycolysis$REACTION[9],referenceData = chemicalData,mCharge = "CHARGE",ids = "NAME"),TRUE)
  checkEquals(isBalanced(reactionList = glycolysis$REACTION[9],referenceData = chemicalData,mWeight = "MASS",ids = "NAME"),TRUE)
}

