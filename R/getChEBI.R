#' @export getChEBI
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Download the ChEBI database
#' @description This function downloads the compounds, formulas, masses and charges from the selected release of the ChEBI database. 
#' The ChEBI database (Chemical Entities of Biological Interest), is a database and ontology of molecular entities focused on 'small' chemical compounds.
#' @param release A release number of the ChEBI database version to be downloaded, by default \code{'latest'} release is downloaded.
#' @return A data.frame with the following data associated to the ChEBI compounds: \itemize{
#' \item \code{"ID"}: The unique identifer
#' \item \code{"ChEBI"}: The name recommended for use in biological databases
#' \item \code{"KEGG"}: The associated name(s) in the KEGG database
#' \item \code{"IUPAC"}: The name(s) generated according to recommendations of IUPAC
#' \item \code{"MetaCyc"}: The associated name(s) in the MetaCyc database
#' \item \code{"ChEMBL"}: The associated name(s) in the ChEMBL database
#' \item \code{"FORMULA"}: The molecular formula
#' \item \code{"MASS"}: The molecular mass
#' \item \code{"MONOISOTOPIC"}: The molecular monoisotopic mass
#' \item \code{"CHARGE"}: The molecular net charge 
#'  }
#' @examples 
#' \dontrun{
#' ChEBI <- getChEBI(
#'  release = "142"
#'  )
#'  }
#' @seealso The ChEBI database webpage: https://www.ebi.ac.uk/chebi/
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
  message("DONE",appendLF = TRUE)
  message("Downloading synonyms ... ",appendLF = FALSE)
  download.file(paste0(ftp,"names.tsv.gz"),paste0(chebi_download,"names.tsv"),quiet = TRUE)
  names <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,"names.tsv"))))
  message("DONE",appendLF = TRUE)
  message("Downloading formulas ... ",appendLF = FALSE)
  download.file(paste0(ftp,"chemical_data.tsv"),paste0(chebi_download,"formulas.tsv"),quiet = TRUE)
  formulas <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,"formulas.tsv"))))
  message("DONE",appendLF = TRUE)
  
  # Building database
  message("Building ChEBI ... ",appendLF = TRUE)
  compoundsParent <- compounds
  compounds[compounds[,"PARENT_ID"]!="null","ID"] <- compounds[compounds[,"PARENT_ID"]!="null","PARENT_ID"]
  compounds <- unique(rbind(compounds,compoundsParent))
  DB <- suppressWarnings((merge(compounds[,c("ID","PARENT_ID","NAME")],names[,c("COMPOUND_ID","SOURCE","NAME")],by.x = "ID",by.y = "COMPOUND_ID",all.x = TRUE)))
  DB[DB[,"PARENT_ID"]!="null","ID"] <- DB[DB[,"PARENT_ID"]!="null","PARENT_ID"]
  DB <- unique(DB)
  #  Associations
  ## CHEBI
  ChEBI <- unique(DB[DB[,"NAME.x"]!="null",c("ID","NAME.x")])
  ChEBI <- ChEBI[complete.cases(ChEBI),]
  colnames(ChEBI) <- c("ID","ChEBI")
  ## KEGG
  message(" KEGG Associations ... ",appendLF = FALSE)
  KEGG <- unique(DB[DB[,"SOURCE"]=="KEGG COMPOUND",c("ID","NAME.y")])
  KEGG <- KEGG[complete.cases(KEGG),]
  colnames(KEGG) <- c("ID","KEGG")
  message("DONE",appendLF = TRUE)
  ## IUPAC
  message(" IUPAC Associations ... ",appendLF = FALSE)
  IUPAC <- unique(DB[DB[,"SOURCE"]=="IUPAC",c("ID","NAME.y")])
  IUPAC <- IUPAC[complete.cases(IUPAC),]
  colnames(IUPAC) <- c("ID","IUPAC")
  message("DONE",appendLF = TRUE)
  ## METACYC
  message(" MetaCyc Associations ... ",appendLF = FALSE)
  MetaCyc <- unique(DB[DB[,"SOURCE"]=="MetaCyc",c("ID","NAME.y")])
  MetaCyc <- MetaCyc[complete.cases(MetaCyc),]
  colnames(MetaCyc) <- c("ID","MetaCyc")
  message("DONE",appendLF = TRUE)
  ## CHEMBL
  message(" ChEMBL Associations ... ",appendLF = FALSE)
  ChEMBL <- unique(DB[DB[,"SOURCE"]=="ChEMBL",c("ID","NAME.y")])
  ChEMBL <- ChEMBL[complete.cases(ChEMBL),]
  colnames(ChEMBL) <- c("ID","ChEMBL")
  message("DONE",appendLF = TRUE)
  DB <- unique(merge(DB["ID"],ChEBI,by = "ID",all.x = TRUE))
  DB <- unique(merge(DB,KEGG,by = "ID",all.x = TRUE))
  DB <- unique(merge(DB,IUPAC,by = "ID",all.x = TRUE))
  DB <- unique(merge(DB,MetaCyc,by = "ID",all.x = TRUE))
  DB <- unique(merge(DB,ChEMBL,by = "ID",all.x = TRUE))
  ## FORMULA
  if("FORMULA" %in% unique(formulas[,"TYPE"])){
    message(" Formula Associations ... ", appendLF = FALSE)
    formula <- formulas[formulas[,"TYPE"]=="FORMULA",c("COMPOUND_ID","CHEMICAL_DATA")]
    colnames(formula) <- c("ID","FORMULA")
    DB <- unique(merge(DB,formula, by = "ID",all.x = TRUE))
    message("DONE",appendLF = TRUE)
    } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  ## MASS
  message("Downloading molecular weights ... ",appendLF = FALSE)
  if("MASS" %in% unique(formulas[,"TYPE"])){
  mass <- formulas[formulas[,"TYPE"]=="MASS",c("COMPOUND_ID","CHEMICAL_DATA")]
  colnames(mass) <- c("ID","MASS")
  DB <- merge(DB,mass, by = "ID",all.x = TRUE)
  message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  message("Downloading monoisotopic molecular weights ... ",appendLF = FALSE)
  ## MMASS
  if("MONOISOTOPIC MASS" %in% unique(formulas[,"TYPE"])){
    mmass <- formulas[formulas[,"TYPE"]=="MONOISOTOPIC MASS",c("COMPOUND_ID","CHEMICAL_DATA")]
    colnames(mmass) <- c("ID","MONOISOTOPIC")
    DB <- merge(DB,mmass, by = "ID",all.x = TRUE)
    message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  message("Downloading molecular charges ... ",appendLF = FALSE)
  ## CHARGE
  if("CHARGE" %in% unique(formulas[,"TYPE"])){
  charge <- formulas[formulas[,"TYPE"]=="CHARGE",c("COMPOUND_ID","CHEMICAL_DATA")]
  colnames(charge) <- c("ID","CHARGE")
  DB <- merge(DB,charge, by = "ID",all.x = TRUE)
  message("DONE",appendLF = TRUE)
  } else { message("NOT AVAILABLE FOR THIS RELEASE")}
  ChEBI <- unique(DB)
  message("DONE",appendLF = TRUE)
  return(ChEBI)
}