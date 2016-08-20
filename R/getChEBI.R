# getChEBI
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

getChEBI <- function(release="latest"){
  # Download folder
  chebi_download <- tempdir()
  # Download releases
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/",paste0(chebi_download,"releases.txt"),quiet = TRUE)
  releases <- gsub("rel","",read.table(paste0(chebi_download,"releases.txt"), quote="\"", comment.char="")[,9])
  message("Validating ChEBI release number ... ",appendLF = FALSE)
  # Match release
  if (release == "latest"){
    release <- max(releases)
  } else {
    release <- releases[match(release,releases)]
  }
  message("OK")
  # Download files
  ftp <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel",release,"/Flat_file_tab_delimited/")
  message("Downloading compounds ... ",appendLF = FALSE)
  download.file(paste0(ftp,"compounds.tsv.gz"),paste0(chebi_download,"compounds.tsv"),quiet = TRUE)
  compounds <- as.data.frame.array(read.delim2(paste0(chebi_download,"compounds.tsv")))
  compounds <- compounds[compounds["STAR"]>=3,c("ID","NAME")]
  compounds[compounds[,"NAME"]=="null","NAME"] <- NA
  message("DONE",appendLF = TRUE)
  message("Downloading synonyms ... ",appendLF = FALSE)
  download.file(paste0(ftp,"names.tsv.gz"),paste0(chebi_download,"names.tsv"),quiet = TRUE)
  names <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,"names.tsv"))))
  kegg <- names[names[,"SOURCE"]=="KEGG COMPOUND",]
  kegg <- names[names[,"TYPE"]=="NAME",c("COMPOUND_ID","NAME")]
  colnames(kegg) <- c("ID","KEGG")
  names <- as.data.frame.array(names[names[,"COMPOUND_ID"]%in%compounds[,"ID"],c("COMPOUND_ID","NAME")])
  colnames(names) <- c("ID","SYNONYMS")
  DB <- merge(compounds,names, by = "ID")
  DB <- merge(DB, kegg, by = "ID")
  message("DONE",appendLF = TRUE)
  message("Downloading formulas ... ",appendLF = FALSE)
  download.file(paste0(ftp,"chemical_data.tsv"),paste0(chebi_download,"formulas.tsv"),quiet = TRUE)
  formulas <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,"formulas.tsv"))))
  if("FORMULA" %in% unique(formulas[,"TYPE"])){
    formula <- formulas[formulas[,"TYPE"]=="FORMULA",c("COMPOUND_ID","CHEMICAL_DATA")]
    colnames(formula) <- c("ID","FORMULA")
    DB <- merge(DB,formula, by = "ID")
    message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  message("Downloading molecular weights ... ",appendLF = FALSE)
  if("MASS" %in% unique(formulas[,"TYPE"])){
  mass <- formulas[formulas[,"TYPE"]=="MASS",c("COMPOUND_ID","CHEMICAL_DATA")]
  colnames(mass) <- c("ID","MASS")
  DB <- merge(DB,mass, by = "ID",all.x = TRUE)
  message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  message("Downloading molecular charges ... ",appendLF = FALSE)
  if("CHARGE" %in% unique(formulas[,"TYPE"])){
  charge <- formulas[formulas[,"TYPE"]=="CHARGE",c("COMPOUND_ID","CHEMICAL_DATA")]
  colnames(charge) <- c("ID","CHARGE")
  DB <- merge(DB,charge, by = "ID",all.x = TRUE)
  message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  message("Downloading monoisotopic molecular weights ... ",appendLF = FALSE)
  if("MONOISOTOPIC MASS" %in% unique(formulas[,"TYPE"])){
  mmass <- formulas[formulas[,"TYPE"]=="MONOISOTOPIC MASS",c("COMPOUND_ID","CHEMICAL_DATA")]
  colnames(mmass) <- c("ID","MONOISOTOPIC")
  DB <- merge(DB,mmass, by = "ID",all.x = TRUE)
  message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  # Building database
  message("Building ChEBI ... ",appendLF = FALSE)
  ChEBI <- unique(DB)
  message("DONE",appendLF = TRUE)
  return(ChEBI)
}