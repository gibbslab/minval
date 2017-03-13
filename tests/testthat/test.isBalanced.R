# Output Class
chemicalData <- read.csv2(system.file("extdata", "chemData.csv", package = "minval"))
test_that("isBalanced function: Output class is wrong", {
  expect_true(is.vector(isBalanced(reactionList = "H2O[c] + H+[c] => ADP[c]",referenceData = chemicalData,ids = "NAME",mFormula = "FORMULA")))
})

# Output Value
test_that("isBalanced function: Output value is wrong", {
  expect_false(isBalanced(reactionList = "H2O[c] + H+[c] => ADP[c]",referenceData = chemicalData,ids = "NAME",mFormula = "FORMULA"))
  expect_false(isBalanced(reactionList = "5 H+[c] => H+[c]",referenceData = chemicalData,ids = "NAME",mFormula = "FORMULA"))
  expect_false(isBalanced(reactionList = "5 H+[c] => H+[c]",referenceData = chemicalData,ids = "NAME",mWeight = "MASS"))
  expect_true(isBalanced(reactionList = "ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]",referenceData = chemicalData,ids = "NAME",mFormula = "FORMULA"))
  expect_true(isBalanced(reactionList = "H+[c] => H+[m]",referenceData = chemicalData,ids = "NAME",mFormula = "FORMULA"))
})