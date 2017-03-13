# Output Class
test_that("isValidSyntax function: Output class is wrong", {
  expect_true(is.vector(isValidSyntax(
    reactionList = c("4A[c] + 3 B[m] => 2.165557 C[e]", "A[c] + B[m] => C[e]")
  )))
})

# Output Value
test_that("isValidSyntax function: Output value is wrong", {
  expect_true(isValidSyntax(reactionList = "A[c] + B[m] => C[e]"))
  expect_true(isValidSyntax(reactionList = "4A[c] + 3 B[m] => 2.165557 C[e]"))
  expect_warning(isValidSyntax(reactionList = "4A[c] + 3 3 B[m] => 2.165557 C[e]"))
  expect_false(suppressWarnings(isValidSyntax(reactionList = "4A[c] + 3 3 B[m] => 2.165557 C[e]")))
})