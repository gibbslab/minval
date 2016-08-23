# convert2SBML
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

convert2SBML<-function(data,outfile){
  # Creating SBMLR model
  model <- convert2SBMLR(data)
  # Writing model
  # This function is a modified copy of saveSBML function included in SBMLR package.
  # Original 'saveSBML' function was writed by Tomas Radivoyevitch
  # Please see: https://www.bioconductor.org/packages/release/bioc/html/SBMLR.html
  writeSBML <- function(model,outfile) {
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
    cat(sprintf("<model id=\"%s\" name=\"%s\">",id,id), file=fid, sep="\n")
    if(nNotes>0){
      cat("<notes>", file=fid, sep="\n")
      cat(" <body xmlns=\"http://www.w3.org/1999/xhtml\">", file=fid, sep="\n")
      sapply(1:nNotes, function(i){cat(sprintf("   <p> %s  </p>",notes[[i]]), file=fid, sep="\n")})
      cat(" </body>", file=fid, sep="\n")
      cat("</notes>", file=fid, sep="\n")
    }
    if(nCompartments>0){
      cat("<listOfCompartments>", file=fid, sep="\n")
      sapply(1:nCompartments, function(i){cat(sprintf("   <compartment id=\"%s\"  name=\"%s\"/>", compartments[[i]][["id"]],compartments[[i]][["name"]]), file=fid, sep="\n")})
      cat("</listOfCompartments>", file=fid, sep="\n")
    }
    if(nSpecies>0){
      cat("<listOfSpecies>", file=fid, sep="\n")
      sapply(1:nSpecies,function(i){cat(sprintf("   <species id=\"%s\"  name=\"%s\"  compartment=\"%s\"/>",species[[i]][["id"]],species[[i]][["name"]],species[[i]][["compartment"]]), file=fid, sep="\n")})
      cat("</listOfSpecies>", file=fid, sep="\n")
    }
    cat("<listOfReactions>", file=fid, sep="\n")
    sapply(1:nReactions, function(i){
      if (is.null(reactions[[i]][["reversible"]])){ print("Internal SBMLR object should have reverse flag set")} else{
        cat(sprintf("  <reaction id=\"%s\"  reversible=\"%s\">",reactions[[i]][["id"]], ifelse(reactions[[i]][["reversible"]],"true","false")), file=fid, sep="\n") 
      }
      gpr = reactions[[i]][["notes"]]
      if(!is.na(gpr[["GPR"]])){
        cat(sprintf("    <notes>"),file=fid,sep="\n")
        cat(sprintf("      <html xmlns=\"http://www.w3.org/1999/xhtml\"><p>GENE_ASSOCIATION: %s</p></html>",gpr[["GPR"]]),file=fid,sep="\n")
        cat(sprintf("    </notes>"),file=fid,sep="\n")
      }
      reactants=reactions[[i]][["reactants"]]
      if (!is.null(reactants[[1]])) {
        cat("    <listOfReactants>", file=fid, sep="\n")
        sapply(1:length(reactants[["reactants"]]),function(j){ cat(sprintf("      <speciesReference species=\"%s\" stoichiometry=\"%s\"/>", reactants[["reactants"]][[j]],reactants[["stoichiometry"]][[j]]), file=fid, sep="\n")})
        cat("    </listOfReactants>", file=fid, sep="\n")
      }
      
      # Just switched the order of these two blocks to fix the errors
      products=reactions[[i]][["products"]]
      if (!is.null(products[[1]])) {
        cat("    <listOfProducts>", file=fid, sep="\n")
        sapply(1:length(products[["products"]]), function(j){ cat(sprintf("      <speciesReference species=\"%s\" stoichiometry=\"%s\"/>", products[["products"]][[j]],products[["stoichiometry"]][[j]]), file=fid, sep="\n")})
        cat("    </listOfProducts>", file=fid, sep="\n")     
      }
      
      cat("  <kineticLaw>", file=fid, sep="\n")
      cat("    <math xmlns=\"http://www.w3.org/1998/Math/MathML\">", file=fid, sep="\n")
      # cat(reactions[[i]]$mathmlLaw)
      ml=XML::saveXML(reactions[[i]]$mathmlLaw,prefix=NULL,file=fid) # annoying warnings were coming from here without file=fid
      cat("    </math>", file=fid, sep="\n")
      
      parameters=reactions[[i]][["parameters"]]
      nlocalParameters = length(parameters)
      if(nlocalParameters > 0){
        cat("    <listOfParameters>", file=fid, sep="\n")
        sapply(1:nlocalParameters,function(j){ cat(sprintf("      <parameter id=\"%s\" value=\"%g\"/>", names(parameters)[j],parameters[[j]]), file=fid, sep="\n")})
        cat("    </listOfParameters>", file=fid, sep="\n")
        cat("    </kineticLaw>", file=fid, sep="\n")
        cat("  </reaction>", file=fid, sep="\n")}})
    
    cat("</listOfReactions>", file=fid, sep="\n")
    cat("</model>", file=fid, sep="\n")
    cat("</sbml>", file=fid, sep="\n")
    close(fid)
  }
  writeSBML(model,outfile)
}