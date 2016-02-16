# minval-internal
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

.metname <- function(met, rm.coef = FALSE) {
  if (rm.coef == FALSE) {
    strsplit(met, paste("[", compartments(met), "]", sep = ""), fixed = TRUE)[[1]]
  }
  else {
    sub("^[0-9]+ ", "", strsplit(met, paste("[", compartments(met), "]", sep = ""), fixed = TRUE)[[1]])
  }
}

.formula2matrix <- function(formula) {
  byatomtype <-
    unlist(regmatches(formula, gregexpr("([A-Z]{1}[a-z]?)([0-9]*)", formula)))
  atomtype <- sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', byatomtype)
  atomnumber <-
    as.numeric(regmatches(byatomtype, gregexpr('[0-9]+', byatomtype)))
  atomnumber[is.na(atomnumber)] <- 1
  tapply(atomnumber, atomtype, sum)
}

.coeficients <- function(met) {
  regmatches(met, gregexpr('^[0-9,.]', met))
}

# .rxnFromModel <- function(file){
#   require(gdata)
#   data <- gdata::read.xls(file , sheet = 1)
#   as.vector(data[-c(data[,1]=="#"),4])
# }

