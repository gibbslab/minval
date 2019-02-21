#' @export downloadChEBI
#' @importFrom stats complete.cases setNames
#' @importFrom utils download.file read.delim2 read.table
#' @author Daniel Camilo Osorio <dcosorioh@tamu.edu>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Download the ChEBI database
#' @description This function downloads the compounds, formulas, masses and charges from the selected release of the ChEBI database.
#' The ChEBI database (Chemical Entities of Biological Interest), is a database and ontology of molecular entities focused on 'small' chemical compounds.
#' @param release A character string with the release number of the ChEBI database version to be downloaded, by default the \code{'latest'} release is downloaded.
#' @param woAssociations A logical value \code{'TRUE'} or \code{'FALSE'} if a light version of the ChEBI database without associations should be returned.
#' @return A data.frame with the following data associated to the ChEBI compounds: \itemize{
#' \item \code{'ID'}: The unique identifer
#' \item \code{'ChEBI'}: The name recommended for use in biological databases
#' \item \code{'KEGG'}: The associated name(s) in the KEGG database
#' \item \code{'IUPAC'}: The name(s) generated according to recommendations of IUPAC
#' \item \code{'MetaCyc'}: The associated name(s) in the MetaCyc database
#' \item \code{'ChEMBL'}: The associated name(s) in the ChEMBL database
#' \item \code{'FORMULA'}: The molecular formula
#' \item \code{'MASS'}: The molecular mass
#' \item \code{'MONOISOTOPIC'}: The molecular monoisotopic mass
#' \item \code{'CHARGE'}: The molecular net charge
#' }
#' If woAssociations is \code{'TRUE'} a A data.frame with the following data is returned: \itemize{
#' \item \code{'NAME'}: The name(s) associated in several biological databases
#' \item \code{'FORMULA'}: The molecular formula
#' \item \code{'MASS'}: The molecular mass
#' \item \code{'MONOISOTOPIC'}: The molecular monoisotopic mass
#' \item \code{'CHARGE'}: The molecular net charge
#' }
#' @examples
#' \dontrun{
#' # Download ChEBI database with associations
#' ChEBI <- downloadChEBI(release = '142')
#'
#' # Download ChEBI database without associations
#' ChEBI <- downloadChEBI(release = '142', woAssociations = TRUE)
#'  }
#' @seealso The ChEBI database webpage: https://www.ebi.ac.uk/chebi/
downloadChEBI <- function(release = "latest",
                          woAssociations = FALSE) {
  # Download folder
  chebi_download <- tempdir()
  # Download releases
  download.file(
    "ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/",
    paste0(chebi_download, "releases.txt"),
    quiet = TRUE, method = "libcurl"
  )
  releases <-
    gsub("rel", "", read.table(
      paste0(chebi_download, "releases.txt"),
      quote = "\"",
      comment.char = ""
    )[, 9])
  message("Validating ChEBI release number ... ", appendLF = FALSE)
  # Match release
  if (release == "latest") {
    release <- max(releases)
  } else {
    release <- releases[match(release, releases)]
  }
  message("OK")
  # Download files
  ftp <-
    paste0(
      "ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel",
      release,
      "/Flat_file_tab_delimited/"
    )
  message("Downloading compounds ... ", appendLF = FALSE)
  download.file(paste0(ftp, "compounds.tsv.gz"),
                paste0(chebi_download, "compounds.tsv"),
                quiet = TRUE, method = "libcurl")
  compounds <- suppressWarnings(
    as.data.frame.array(read.delim2(paste0(chebi_download, "compounds.tsv"))))
  message("DONE", appendLF = TRUE)
  message("Downloading synonyms ... ", appendLF = FALSE)
  download.file(paste0(ftp, "names.tsv.gz"),
                paste0(chebi_download, "names.tsv"),
                quiet = TRUE, method = "libcurl")
  names <-
    suppressWarnings(as.data.frame.array(read.delim2(paste0(
      chebi_download, "names.tsv"
    ))))
  message("DONE", appendLF = TRUE)
  message("Downloading formulas ... ", appendLF = FALSE)
  download.file(paste0(ftp, "chemical_data.tsv"),
                paste0(chebi_download, "formulas.tsv"),
                quiet = TRUE, method = "libcurl")
  formulas <-
    suppressWarnings(as.data.frame.array(read.delim2(
      paste0(chebi_download, "formulas.tsv")
    )))
  message("DONE", appendLF = TRUE)
  
  # Building database
  message("Building ChEBI ... ", appendLF = TRUE)
  compounds <- compounds[compounds[, "STAR"] >= 3, ]
  latest <- compounds[, c("ID", "NAME")]
  old <- compounds[, c("ID", "PARENT_ID")]
  old <- merge(
    x = old,
    y = latest,
    by.x = "PARENT_ID",
    by.y = "ID"
  )
  compounds <- rbind(latest, old[, c("ID", "NAME")])
  
  if(any(sapply(compounds[, "NAME"] == "null", isTRUE))){
    compounds[compounds[, "NAME"] == "null", "NAME"] <- NA
  }
  compounds <- compounds[complete.cases(compounds), ]
  DB <-
    suppressWarnings((
      merge(
        compounds[, c("ID", "NAME")],
        names[, c("COMPOUND_ID", "SOURCE", "NAME")],
        by.x = "ID",
        by.y = "COMPOUND_ID",
        all.x = TRUE
      )
    ))
  #  Associations
  ## CHEBI
  ChEBI <- unique(DB[, c("ID", "NAME.x")])
  colnames(ChEBI) <- c("ID", "ChEBI")
  ## KEGG
  message(" KEGG Associations ... ", appendLF = FALSE)
  KEGG <-
    unique(DB[DB[, "SOURCE"] == "KEGG COMPOUND", c("ID", "NAME.y")])
  KEGG <- KEGG[complete.cases(KEGG), ]
  colnames(KEGG) <- c("ID", "KEGG")
  message("DONE", appendLF = TRUE)
  ## IUPAC
  message(" IUPAC Associations ... ", appendLF = FALSE)
  IUPAC <- unique(DB[DB[, "SOURCE"] == "IUPAC", c("ID", "NAME.y")])
  IUPAC <- IUPAC[complete.cases(IUPAC), ]
  colnames(IUPAC) <- c("ID", "IUPAC")
  message("DONE", appendLF = TRUE)
  ## METACYC
  message(" MetaCyc Associations ... ", appendLF = FALSE)
  MetaCyc <- unique(DB[DB[, "SOURCE"] == "MetaCyc", c("ID", "NAME.y")])
  MetaCyc <- MetaCyc[complete.cases(MetaCyc), ]
  colnames(MetaCyc) <- c("ID", "MetaCyc")
  message("DONE", appendLF = TRUE)
  ## CHEMBL
  message(" ChEMBL Associations ... ", appendLF = FALSE)
  ChEMBL <- unique(DB[DB[, "SOURCE"] == "ChEMBL", c("ID", "NAME.y")])
  ChEMBL <- ChEMBL[complete.cases(ChEMBL), ]
  colnames(ChEMBL) <- c("ID", "ChEMBL")
  message("DONE", appendLF = TRUE)
  DB <- unique(merge(DB["ID"], ChEBI, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, KEGG, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, IUPAC, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, MetaCyc, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, ChEMBL, by = "ID", all.x = TRUE))
  rm(ChEBI,
     ChEMBL,
     compounds,
     IUPAC,
     KEGG,
     latest,
     MetaCyc,
     names,
     old)
  ## FORMULA
  if ("FORMULA" %in% unique(formulas[, "TYPE"])) {
    message(" Formula Associations ... ", appendLF = FALSE)
    formula <-
      formulas[formulas[, "TYPE"] == "FORMULA", c("COMPOUND_ID", "CHEMICAL_DATA")]
    colnames(formula) <- c("ID", "FORMULA")
    DB <- merge(DB, formula, by = "ID", all.x = TRUE)
    DB <-
      merge(DB, DB[, c("ChEBI", "FORMULA")], by = "ChEBI", all.x = TRUE)
    DB[is.na(DB[, "FORMULA.x"]), "FORMULA.x"] <- "null"
    DB[is.na(DB[, "FORMULA.y"]), "FORMULA.y"] <- "null"
    DB[DB[, "FORMULA.x"] != "null" &
         DB[, "FORMULA.y"] == "null", "FORMULA.y"] <-
      DB[DB[, "FORMULA.x"] != "null" &
           DB[, "FORMULA.y"] == "null", "FORMULA.x"]
    DB[DB[, "FORMULA.y"] != "null" &
         DB[, "FORMULA.x"] == "null", "FORMULA.x"] <-
      DB[DB[, "FORMULA.y"] != "null" &
           DB[, "FORMULA.x"] == "null", "FORMULA.y"]
    DB <-
      unique(DB[DB[, "FORMULA.x"] != "null" &
                  DB[, "FORMULA.y"] != "null", c("ID",
                                                 "ChEBI",
                                                 "KEGG",
                                                 "IUPAC",
                                                 "MetaCyc",
                                                 "ChEMBL",
                                                 "FORMULA.x")])
    rm(formula)
    message("DONE", appendLF = TRUE)
  } else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  ## MASS
  message("Downloading molecular weights ... ", appendLF = FALSE)
  if ("MASS" %in% unique(formulas[, "TYPE"])) {
    mass <-
      formulas[formulas[, "TYPE"] == "MASS", c("COMPOUND_ID", "CHEMICAL_DATA")]
    colnames(mass) <- c("ID", "MASS")
    DB <- merge(DB, mass, by = "ID", all.x = TRUE)
    DB <- merge(DB, DB[, c("ChEBI", "MASS")], by = "ChEBI", all.x = TRUE)
    DB[is.na(DB[, "MASS.x"]), "MASS.x"] <- "null"
    DB[is.na(DB[, "MASS.y"]), "MASS.y"] <- "null"
    DB[DB[, "MASS.x"] != "null" &
         DB[, "MASS.y"] == "null", "MASS.y"] <-
      DB[DB[, "MASS.x"] != "null" & DB[, "MASS.y"] == "null", "MASS.x"]
    DB[DB[, "MASS.y"] != "null" &
         DB[, "MASS.x"] == "null", "MASS.x"] <-
      DB[DB[, "MASS.y"] != "null" & DB[, "MASS.x"] == "null", "MASS.y"]
    DB <-
      unique(DB[, c("ID",
                    "ChEBI",
                    "KEGG",
                    "IUPAC",
                    "MetaCyc",
                    "ChEMBL",
                    "FORMULA.x",
                    "MASS.x")])
    rm(mass)
    message("DONE", appendLF = TRUE)
  } else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  message("Downloading monoisotopic molecular weights ... ", appendLF = FALSE)
  ## MMASS
  if ("MONOISOTOPIC MASS" %in% unique(formulas[, "TYPE"])) {
    mmass <-
      formulas[formulas[, "TYPE"] == "MONOISOTOPIC MASS", c("COMPOUND_ID", "CHEMICAL_DATA")]
    colnames(mmass) <- c("ID", "MONOISOTOPIC")
    DB <- merge(DB, mmass, by = "ID", all.x = TRUE)
    DB <-
      merge(DB, DB[, c("ChEBI", "MONOISOTOPIC")], by = "ChEBI", all.x = TRUE)
    DB[is.na(DB[, "MONOISOTOPIC.x"]), "MONOISOTOPIC.x"] <- "null"
    DB[is.na(DB[, "MONOISOTOPIC.y"]), "MONOISOTOPIC.y"] <- "null"
    DB[DB[, "MONOISOTOPIC.x"] != "null" &
         DB[, "MONOISOTOPIC.y"] == "null", "MONOISOTOPIC.y"] <-
      DB[DB[, "MONOISOTOPIC.x"] != "null" &
           DB[, "MONOISOTOPIC.y"] == "null", "MONOISOTOPIC.x"]
    DB[DB[, "MONOISOTOPIC.y"] != "null" &
         DB[, "MONOISOTOPIC.x"] == "null", "MONOISOTOPIC.x"] <-
      DB[DB[, "MONOISOTOPIC.y"] != "null" &
           DB[, "MONOISOTOPIC.x"] == "null", "MONOISOTOPIC.y"]
    DB <-
      unique(DB[, c(
        "ID",
        "ChEBI",
        "KEGG",
        "IUPAC",
        "MetaCyc",
        "ChEMBL",
        "FORMULA.x",
        "MASS.x",
        "MONOISOTOPIC.x"
      )])
    rm(mmass)
    message("DONE", appendLF = TRUE)
  } else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  message("Downloading molecular charges ... ", appendLF = FALSE)
  ## CHARGE
  if ("CHARGE" %in% unique(formulas[, "TYPE"])) {
    charge <-
      formulas[formulas[, "TYPE"] == "CHARGE", c("COMPOUND_ID", "CHEMICAL_DATA")]
    colnames(charge) <- c("ID", "CHARGE")
    DB <- merge(DB, charge, by = "ID", all.x = TRUE)
    DB <- merge(DB, DB[, c("ChEBI", "CHARGE")], by = "ChEBI", all.x = TRUE)
    DB[is.na(DB[, "CHARGE.x"]), "CHARGE.x"] <- "null"
    DB[is.na(DB[, "CHARGE.y"]), "CHARGE.y"] <- "null"
    DB[DB[, "CHARGE.x"] != "null" &
         DB[, "CHARGE.y"] == "null", "CHARGE.y"] <-
      DB[DB[, "CHARGE.x"] != "null" & DB[, "CHARGE.y"] == "null", "CHARGE.x"]
    DB[DB[, "CHARGE.y"] != "null" &
         DB[, "CHARGE.x"] == "null", "CHARGE.x"] <-
      DB[DB[, "CHARGE.y"] != "null" & DB[, "CHARGE.x"] == "null", "CHARGE.y"]
    DB <-
      unique(DB[, c(
        "ID",
        "ChEBI",
        "KEGG",
        "IUPAC",
        "MetaCyc",
        "ChEMBL",
        "FORMULA.x",
        "MASS.x",
        "MONOISOTOPIC.x",
        "CHARGE.x"
      )])
    message("DONE", appendLF = TRUE)
  } else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  DB[DB == "null"] <- NA
  DB <-
    unique(DB[complete.cases(DB[, c("ID",
                                    "ChEBI",
                                    "FORMULA.x",
                                    "MASS.x",
                                    "MONOISOTOPIC.x",
                                    "CHARGE.x")]), ])
  colnames(DB) <-
    c(
      "ID",
      "ChEBI",
      "KEGG",
      "IUPAC",
      "MetaCyc",
      "ChEMBL",
      "FORMULA",
      "MASS",
      "MONOISOTOPIC",
      "CHARGE"
    )
  if (woAssociations == TRUE) {
    compounds <-
      unique(rbind(
        setNames(DB[, c("ChEBI", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE")], c(
          "NAME", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE"
        )),
        setNames(DB[, c("KEGG", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE")], c(
          "NAME", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE"
        )),
        setNames(DB[, c("IUPAC", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE")], c(
          "NAME", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE"
        )),
        setNames(DB[, c("MetaCyc", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE")], c(
          "NAME", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE"
        )),
        setNames(DB[, c("ChEMBL", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE")], c(
          "NAME", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE"
        ))
      ))
    compounds <- compounds[complete.cases(compounds), ]
    return(compounds)
  } else {
    return(DB)
  }
}
