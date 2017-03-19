writeTSV <-
  function (modelData,
            modelID = "model",
            outputFile,
            boundary = "b") {
    if (class(modelData) == "data.frame") {
      # Check valid structure, column names and valid ID's
      modelData <- validateData(modelData = modelData)
      # Remove comments
      modelData <- removeComments(modelData = modelData)
      # Validate stoichiometric syntax
      modelData <- modelData[validateSyntax(modelData[["REACTION"]]), ]
    } else if (class(modelData) == "modelorg") {
      modelData <- convertData(model = modelData)
    } else {
      stop("Input format not supported.")
    }
    # writeTSV _react.tsv
    outputData <- NULL
    outputData$abbreviation <- modelData[["ID"]]
    outputData$name <- modelData[["DESCRIPTION"]]
    outputData$equation <-
      rearmReactions(
        S = stoichiometricMatrix(modelData[["REACTION"]]),
        reversible = grepl("<=>", modelData[["REACTION"]]),
        type = "TSV",
        boundary = boundary
      )
    outputData$reversible <-
      ifelse(test = grepl("<=>", modelData[["REACTION"]]),
             yes = "reversible",
             no = "irreversible")
    outputData$compartment <-
      unlist(lapply(modelData[["REACTION"]], function(reaction) {
        paste0(compartments(reaction), collapse = " , ")
      }))
    outputData$lowbnd <- as.character(modelData[["LOWER.BOUND"]])
    outputData$uppbnd <- as.character(modelData[["UPPER.BOUND"]])
    outputData$obj_coef <- as.character(modelData[["OBJECTIVE"]])
    outputData$rule <- modelData[["GPR"]]
    outputData$subsystem <- rep("", length(modelData[["REACTION"]]))
    outputData <- data.frame(outputData)
    write.table(
      outputData,
      file = paste0(outputFile, "_react.tsv"),
      quote = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
    
    # writeTSV _met.tsv
    outputData <- NULL
    outputData$abbreviation <-
      gsub("\\+", "_ChargedP", gsub(
        "[[:space:]]+",
        "_",
        metabolites(modelData[["REACTION"]], woCompartment = TRUE)
      ))
    outputData$name <-
      metabolites(modelData[["REACTION"]], woCompartment = TRUE)
    outputData$compartment <-
      compartments(metabolites(modelData[["REACTION"]]), uniques = FALSE)
    outputData <- data.frame(outputData)
    write.table(
      outputData,
      file = paste0(outputFile, "_met.tsv"),
      quote = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
    
    # writeTSV  _desc.tsv
    outputData <- NULL
    outputData$name <- modelID
    outputData$id <- paste(modelID, date(), sep = " - ")
    write.table(
      outputData,
      file = paste0(outputFile, "_desc.tsv"),
      quote = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }
