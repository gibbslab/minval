# minval-internal
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

.remove.compartment <- function(met, rm.coef = FALSE) {
  met <- gsub("^[[:blank:]]*","",met)
  met <- gsub("[[:blank:]]*$","",met)
  if (rm.coef == TRUE) {
    met <- gsub("^[[:digit:]][[:graph:]]*[[:blank:]]","",met)
  }
  gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$","",met)
}

.formula2matrix <- function(formula) {
  byatomtype <- unlist(regmatches(formula, gregexpr("([[:alpha:]]{1}[[:alpha:]]*)([0-9]*)", formula)))
  atomtype <- sub("([A-Z]{1}[a-z]?)([0-9]*)", '\\1', byatomtype)
  atomnumber <- as.numeric(regmatches(byatomtype, gregexpr('[0-9]+', byatomtype)))
  atomnumber[is.na(atomnumber)] <- 1
  tapply(atomnumber, atomtype, sum)
}

.remove.coefficients <- function(met){
  met <- gsub("^[[:digit:]]+[[:punct:]]*[[:digit:]]*[[:blank:]]+","",met)
  return(met)
}

.coefficients <- function(met) {
  met <- regmatches(met, gregexpr('^[[:digit:]][[:punct:]]*[[:digit:]]*[[:blank:]]+', met))
  met[lengths(met)==0] <- 1
  met <- gsub("[[:blank:]]*$","",met)
  met <- as.numeric(met)
  return(met)
}

.atoms <- function(metabolites) {
  coef <- as.numeric(sapply(metabolites, .coefficients))
  formula <- metabolites(metabolites)
  unlist(mapply(function(coef, formula){rep(formula,coef)}, coef = coef, formula = formula,SIMPLIFY = FALSE))
}

.remove.spaces <- function(metabolite){
  metabolite <- gsub("^[[:space:]]","",metabolite)
  metabolite <- gsub("[[:space:]]$","",metabolite)
  return(metabolite)
}

.get.right <- function(reaction){
  .remove.spaces(unlist(strsplit(unlist(strsplit(reaction,"[[:blank:]]*<?=>[[:blank:]]*"))[2],"[[:blank:]]+\\+[[:blank:]]+")))
}

.get.left <- function(reaction){
  .remove.spaces(unlist(strsplit(unlist(strsplit(reaction,"[[:blank:]]*<?=>[[:blank:]]*"))[1],"[[:blank:]]+\\+[[:blank:]]+")))
}

.join.cm <- function(coefficient,metabolite){
  joined <- paste(mapply(function(coefficient,metabolite){paste(coefficient,metabolite,collapse =" ")},coefficient=coefficient,metabolite=metabolite),collapse = " + ")
  return(joined)
}

.join.reaction <- function(reactant,reversibility,product){
  if (reversibility) {
    reaction <- paste(reactant,product,sep = " <=> ")
  } else {
    reaction <- paste(reactant,product,sep = " => ")
  }
}

.validate.xls <- function(data){
  names <- colnames(data)
  if (length(grep("^ID$",names,ignore.case = TRUE))==0){stop("Reaction ID's not found")}
  if (!identical(data[,"ID"],unique(data[,"ID"]))){stop("Reaction ID's must be unique")}
  if (length(grep("^EQUATION$",names,ignore.case = TRUE))==0){stop("Equations not found")}
  if (length(grep("^LOWER.BOUND$",names,ignore.case = TRUE))==0){stop("LB not found")}
  if (length(grep("^UPPER.BOUND$",names,ignore.case = TRUE))==0){stop("UB not found")}
  if (length(grep("^OBJECTIVE$",names,ignore.case = TRUE))==0){stop("OBJECTIVE not found")}
}

.remove.comments <- function(data){
  if (length(grep("#",data[,1]))>0){
    data <- data[!data[,1]=="#"]
  }
  return (data)
}

.fill.reactions <- function(rxnid,data){
  LB = ifelse(is.na(data[data[,"ID"]%in%rxnid,"LOWER.BOUND"]),ifelse(grepl("<=>",data[data[,"ID"]%in%rxnid,"EQUATION"]),-1000,0),data[data[,"ID"]%in%rxnid,"LOWER.BOUND"])
  UB = ifelse(is.na(data[data[,"ID"]%in%rxnid,"UPPER.BOUND"]),1000,data[data[,"ID"]%in%rxnid,"UPPER.BOUND"])
  rev <- grepl("<=>",data[data[,"ID"]%in%rxnid,"EQUATION"])
  left <- .get.left(data[data[,"ID"]%in%rxnid,"EQUATION"])
  right <- .get.right(data[data[,"ID"]%in%rxnid,"EQUATION"])
  reac<-list(id=as.vector(rxnid), 
             reversible=rev,
             reactants=list(reactants=metabolites(left),stoichiometry=.coefficients(left)),
             products=list(products=ifelse(is.na(right),paste0(metabolites(left,woCompartment=TRUE),"[b]"),metabolites(right)),stoichiometry=.coefficients(right)),
             parameters=c(LOWER_BOUND = LB,
                          UPPER_BOUND = UB,
                          OBJECTIVE_COEFFICIENT = ifelse(is.na(data[data[,"ID"]%in%rxnid,"OBJECTIVE"]),0,data[data[,"ID"]%in%rxnid,"OBJECTIVE"]),
                          FLUX_VALUE = 0),
             mathmlLaw = xmlNode("ci","FLUX_VALUE"),
             strlaw = "FLUX_VALUE")
}


# write.xml
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# This function is a modified copy of saveSBML function included in SBMLR package.
# Original 'saveSBML' function was writed by Tomas Radivoyevitch
# Please see: https://www.bioconductor.org/packages/release/bioc/html/SBMLR.html

.write.xml <- function(model,outfile) {
  fid <- file(outfile, "w")
  sbml=model[["sbml"]]
  id=model[["id"]]
  notes=model[["notes"]]
  htmlNotes=model[["htmlNotes"]]
  compartments=model[["compartments"]]
  species=model[["species"]]
  globalParameters=model[["globalParameters"]]
  rules=model[["rules"]]
  reactions=model[["reactions"]]
  units=model[["units"]]
  
  nNotes=length(notes)
  nCompartments=length(compartments)
  nReactions=length(reactions)
  nSpecies=length(species)
  nGlobalParameters=length(globalParameters)
  nRules=length(rules)
  nUnits=length(units)
  
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file=fid, sep="\n")
  cat("<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">", file=fid, sep="\n")
  cat(sprintf("<model id=\"%s\">",id), file=fid, sep="\n")
  if(nNotes>0){
    cat("<notes>", file=fid, sep="\n")
    cat(" <body xmlns=\"http://www.w3.org/1999/xhtml\">", file=fid, sep="\n")
    sapply(seq_along(nNotes), function(i){cat(sprintf("   <p> %s  </p>",notes[[i]]), file=fid, sep="\n")})
    cat(" </body>", file=fid, sep="\n")
    cat("</notes>", file=fid, sep="\n")
  }
  if(nCompartments>0){
    cat("<listOfCompartments>", file=fid, sep="\n")
    sapply(seq_along(nCompartments), function(i){cat(sprintf("   <compartment id=\"%s\"  name=\"%s\"/>", compartments[[i]][["id"]],compartments[[i]][["name"]]), file=fid, sep="\n")})
    cat("</listOfCompartments>", file=fid, sep="\n")
  }
  if(nSpecies>0){
    cat("<listOfSpecies>", file=fid, sep="\n")
    sapply(seq_along(nSpecies),function(i){cat(sprintf("   <species id=\"%s\"  name=\"%s\"  compartment=\"%s\"/>",species[[i]][["id"]],species[[i]][["name"]],species[[i]][["compartment"]]), file=fid, sep="\n")})
           cat("</listOfSpecies>", file=fid, sep="\n")
  }
  cat("<listOfReactions>", file=fid, sep="\n")
  sapply(seq_along(nReactions), function(i){
    if (is.null(reactions[[i]][["reversible"]]))
      print("Internal SBMLR object should have reverse flag set") else
        cat(sprintf("  <reaction id=\"%s\"  reversible=\"%s\">",reactions[[i]][["id"]],
                    ifelse(reactions[[i]][["reversible"]],"true","false")), file=fid, sep="\n")
    reactants=reactions[[i]][["reactants"]]
    if (!is.null(reactants[[1]])) {
      cat("    <listOfReactants>", file=fid, sep="\n")
      sapply(seq_along(length(reactants[["reactants"]])),function(i){ cat(sprintf("      <speciesReference species=\"%s\" stoichiometry=\"%s\"/>", reactants[["reactants"]][[j]],reactants[["stoichiometry"]][[j]]), file=fid, sep="\n")})
      cat("    </listOfReactants>", file=fid, sep="\n")
    }
    
    # Just switched the order of these two blocks to fix the errors
    products=reactions[[i]][["products"]]
    if (!is.null(products[[1]])) {
      cat("    <listOfProducts>", file=fid, sep="\n")
      sapply(seq_along(length(products[["products"]])), function(j){ cat(sprintf("      <speciesReference species=\"%s\" stoichiometry=\"%s\"/>", products[["products"]][[j]],products[["stoichiometry"]][[j]]), file=fid, sep="\n")})
      cat("    </listOfProducts>", file=fid, sep="\n")     
    }
    
    cat("  <kineticLaw>", file=fid, sep="\n")
    cat("    <math xmlns=\"http://www.w3.org/1998/Math/MathML\">", file=fid, sep="\n")
    # cat(reactions[[i]]$mathmlLaw)
    ml=saveXML(reactions[[i]]$mathmlLaw,prefix=NULL,file=fid) # annoying warnings were coming from here without file=fid
    cat("    </math>", file=fid, sep="\n")
    
    parameters=reactions[[i]][["parameters"]]
    nlocalParameters = length(parameters)
    if(nlocalParameters > 0)			#if local parameters exist, we write it. Else we avoid it.
    { cat("    <listOfParameters>", file=fid, sep="\n")
      for (j in 1:nlocalParameters)		# Write each and every local parameter.
        cat(sprintf("      <parameter id=\"%s\" value=\"%g\"/>", names(parameters)[j],parameters[[j]]), file=fid, sep="\n")
      cat("    </listOfParameters>", file=fid, sep="\n")   }
    
    cat("    </kineticLaw>", file=fid, sep="\n")
    cat("  </reaction>", file=fid, sep="\n")})
  cat("</listOfReactions>", file=fid, sep="\n")
  cat("</model>", file=fid, sep="\n")
  cat("</sbml>", file=fid, sep="\n")
  close(fid)
}

.create.model <- function(){
  model <- list(
    id = "",
    notes = c(""),
    compartments = list(),
    species = list(),
    reactions = list(),
    globalParameters = list(),
    rules = list()
  )
  model <-structure(model,class ="SBMLR")
  return(model)
}