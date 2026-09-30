#' Import spectral flow cytometry CSV data
#'
#' Imports spectral flow cytometry CSV files for use with AlgaeScan.
#'
#' Data can be provided in two ways:
#'
#' 1. As separate CSV files containing algae and non-algae events.
#' 2. As a single CSV file containing data that have already been combined.
#'
#' When separate files are provided, AlgaeScan creates a `type` column
#' containing the standardized labels `"algae"` and `"non-algae"`.
#'
#' When a combined file is provided, the user can optionally specify the
#' column containing algae/non-algae labels using `type_col`. These labels
#' are then converted to the standardized AlgaeScan labels.
#'
#' A combined file without known labels can also be imported, for example
#' when analysing new or unknown samples.
#'
#' @param algae_files Character vector containing paths to CSV files
#'   containing algae events. Default is `NULL`.
#'
#' @param non_algae_files Character vector containing paths to CSV files
#'   containing non-algae events. Default is `NULL`.
#'
#' @param combined_file Character string containing the path to a single
#'   CSV file in which the data have already been combined. This cannot
#'   be used together with `algae_files` or `non_algae_files`.
#'   Default is `NULL`.
#'
#' @param type_col Optional character string giving the name of the column
#'   containing algae/non-algae labels when using `combined_file`.
#'   Default is `NULL`.
#'
#' @param algae_label Value in `type_col` corresponding to algae events.
#'   Default is `"algae"`.
#'
#' @param non_algae_label Value in `type_col` corresponding to non-algae
#'   events. Default is `"non-algae"`.
#'
#' @return A `data.table` containing the imported spectral flow cytometry
#'   data. For separate files, `measure` and `type` columns are added.
#'   For combined files, `measure` is added if not already present, and
#'   `type` is created when `type_col` is specified.
#'
#' @examples
#' \dontrun{
#'
#' # Separate algae and non-algae files
#' data <- AlgaeScan_import(
#'   algae_files = c("algae_1.csv", "algae_2.csv"),
#'   non_algae_files = c("background_1.csv", "background_2.csv")
#' )
#'
#' # Already combined and labelled dataset
#' data <- AlgaeScan_import(
#'   combined_file = "combined_data.csv",
#'   type_col = "group",
#'   algae_label = "phytoplankton",
#'   non_algae_label = "background"
#' )
#'
#' # Unlabelled dataset
#' data <- AlgaeScan_import(
#'   combined_file = "unknown_sample.csv"
#' )
#'
#' }
#'
#' @import data.table
#' @export
AlgaeScan_import <- function(
    algae_files = NULL,
    non_algae_files = NULL,
    combined_file = NULL,
    type_col = NULL,
    algae_label = "algae",
    non_algae_label = "non-algae"
) {
  
  # ------------------------------------------------------------
  # Determine import mode
  # ------------------------------------------------------------
  
  separate_mode <- !is.null(algae_files) | !is.null(non_algae_files)
  
  combined_mode <- !is.null(combined_file)
  
  
  # ------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------
  
  if (!separate_mode & !combined_mode) {
    
    stop(
      "Please provide `algae_files` and/or `non_algae_files`, or `combined_file`."
    )
  }
  
  
  if (separate_mode & combined_mode) {
    
    stop(
      "Use either separate algae/non-algae files or `combined_file`, not both."
    )
  }
  
  
  # ============================================================
  # MODE 1: separate algae and non-algae files
  # ============================================================
  
  if (separate_mode) {
    
    all_files <- c(algae_files, non_algae_files)
    
    missing_files <- all_files[!file.exists(all_files)]
    
    if (length(missing_files) > 0) {
      
      stop(
        paste(
          "The following files were not found:",
          paste(missing_files, collapse = "\n"),
          sep = "\n"
        )
      )
    }
    
    
    list_data <- list()
    
    
    # ----------------------------------------------------------
    # Import algae files
    # ----------------------------------------------------------
    
    if (!is.null(algae_files)) {
      
      for (file in algae_files) {
        
        message(
          "Reading algae file: ",
          basename(file)
        )
        
        data <- fread(file)
        
        
        # Keep source file information
        if (!"measure" %in% colnames(data)) {
          
          data[, measure := basename(file)]
        }
        
        
        # Add algae label
        data[, type := "algae"]
        
        
        list_data[[length(list_data) + 1]] <- data
      }
    }
    
    
    # ----------------------------------------------------------
    # Import non-algae files
    # ----------------------------------------------------------
    
    if (!is.null(non_algae_files)) {
      
      for (file in non_algae_files) {
        
        message(
          "Reading non-algae file: ",
          basename(file)
        )
        
        data <- fread(file)
        
        
        # Keep source file information
        if (!"measure" %in% colnames(data)) {
          
          data[, measure := basename(file)]
        }
        
        
        # Add non-algae label
        data[, type := "non-algae"]
        
        
        list_data[[length(list_data) + 1]] <- data
      }
    }
    
    
    # ----------------------------------------------------------
    # Combine files
    # ----------------------------------------------------------
    
    data <- rbindlist(
      list_data,
      use.names = TRUE
    )
  }
  
  
  # ============================================================
  # MODE 2: already combined CSV file
  # ============================================================
  
  if (combined_mode) {
    
    if (length(combined_file) != 1) {
      
      stop(
        "`combined_file` must contain the path to one CSV file."
      )
    }
    
    
    if (!file.exists(combined_file)) {
      
      stop(
        paste(
          "File not found:",
          combined_file
        )
      )
    }
    
    
    message(
      "Reading combined file: ",
      basename(combined_file)
    )
    
    
    data <- fread(combined_file)
    
    
    # ----------------------------------------------------------
    # Add source information if needed
    # ----------------------------------------------------------
    
    if (!"measure" %in% colnames(data)) {
      
      data[, measure := basename(combined_file)]
    }
    
    
    # ----------------------------------------------------------
    # Standardize algae/non-algae labels
    # ----------------------------------------------------------
    
    if (!is.null(type_col)) {
      
      if (!type_col %in% colnames(data)) {
        
        stop(
          paste0(
            "Column `",
            type_col,
            "` was not found in the dataset."
          )
        )
      }
      
      
      original_type <- data[[type_col]]
      
      
      # Check labels
      valid_labels <- original_type %in%
        c(algae_label, non_algae_label)
      
      
      if (any(!valid_labels)) {
        
        stop(
          paste0(
            "Unexpected labels found in `",
            type_col,
            "`. Expected only `",
            algae_label,
            "` and `",
            non_algae_label,
            "`."
          )
        )
      }
      
      
      # Create standardized AlgaeScan type column
      data[
        original_type == algae_label,
        type := "algae"
      ]
      
      data[
        original_type == non_algae_label,
        type := "non-algae"
      ]
    }
  }
  
  
  # ------------------------------------------------------------
  # Import summary
  # ------------------------------------------------------------
  
  message("")
  message("AlgaeScan import complete.")
  message(
    "Events: ",
    format(nrow(data), big.mark = ",")
  )
  
  
  if ("type" %in% colnames(data)) {
    
    message(
      "Algae events: ",
      format(
        sum(data$type == "algae"),
        big.mark = ","
      )
    )
    
    message(
      "Non-algae events: ",
      format(
        sum(data$type == "non-algae"),
        big.mark = ","
      )
    )
  }
  
  
  return(data)
}



#' Load an AlgaeScan model
#'
#' Loads a model stored as an RDS file.
#'
#' Native `AlgaeScan_model` objects are restored directly and can
#' immediately be used with [AlgaeScan_predict()].
#'
#' Legacy `caret::train` objects are also supported. These models are
#' wrapped in an `AlgaeScan_model` object and the available predictor
#' names, class labels, tuning parameters and validation results are
#' recovered when possible.
#'
#' The loader is model-independent and does not require separate logic
#' for Random Forest, KNN, neural networks or other caret models.
#'
#' @param file Character string giving the path of the model RDS file.
#'
#' @param target Optional character string giving the response variable
#'   for a legacy caret model, for example `"type"` or `"Species"`.
#'   Native AlgaeScan models already contain this information.
#'
#' @param label_map Optional named character vector mapping internal
#'   model class labels to desired output labels.
#'
#'   For example:
#'
#'   \code{
#'   c(
#'     "algae" = "algae",
#'     "non_algae" = "non-algae"
#'   )
#'   }
#'
#'   If `NULL`, legacy class labels are preserved exactly as stored
#'   in the caret model.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @examples
#' \dontrun{
#'
#' # Native AlgaeScan model
#' model_A <- AlgaeScan_load_model(
#'   "model_A.rds"
#' )
#'
#' # Legacy caret model
#' model_A_old <- AlgaeScan_load_model(
#'   "old_model_A.rds",
#'   target = "type",
#'   label_map = c(
#'     "algae" = "algae",
#'     "non_algae" = "non-algae"
#'   )
#' )
#'
#' }
#'
#' @export
AlgaeScan_load_model <- function(
    file,
    target = NULL,
    label_map = NULL
) {
  
  # ============================================================
  # 1. Check input file
  # ============================================================
  
  if (
    !is.character(file) ||
    length(file) != 1 ||
    is.na(file) ||
    file == ""
  ) {
    
    stop(
      "`file` must be a single valid file path."
    )
  }
  
  
  file <- path.expand(
    file
  )
  
  
  if (!file.exists(file)) {
    
    stop(
      paste0(
        "Model file does not exist:\n",
        file
      )
    )
  }
  
  
  # ============================================================
  # 2. Check optional arguments
  # ============================================================
  
  if (!is.null(target)) {
    
    if (
      !is.character(target) ||
      length(target) != 1 ||
      is.na(target) ||
      target == ""
    ) {
      
      stop(
        "`target` must be a single character string or NULL."
      )
    }
  }
  
  
  if (!is.null(label_map)) {
    
    if (
      !is.character(label_map) ||
      is.null(names(label_map)) ||
      any(names(label_map) == "")
    ) {
      
      stop(
        "`label_map` must be a named character vector."
      )
    }
  }
  
  
  # ============================================================
  # 3. Read model
  # ============================================================
  
  message(
    "Loading model..."
  )
  
  
  model_object <- tryCatch(
    
    readRDS(
      file
    ),
    
    error = function(e) {
      
      stop(
        paste0(
          "The model could not be loaded.\n",
          "Original error: ",
          conditionMessage(e)
        )
      )
    }
  )
  
  
  # ============================================================
  # 4. Native AlgaeScan model
  # ============================================================
  
  if (
    inherits(
      model_object,
      "AlgaeScan_model"
    )
  ) {
    
    # Native AlgaeScan models already contain the fitted model,
    # predictor names, class labels, training settings and
    # reproducibility information.
    
    if (!is.null(target)) {
      
      model_object$target <- target
    }
    
    
    if (!is.null(label_map)) {
      
      # Determine the internal class names already stored in the
      # AlgaeScan model.
      
      internal_classes <- NULL
      
      
      if (!is.null(model_object$label_map)) {
        
        internal_classes <- names(
          model_object$label_map
        )
      }
      
      
      if (
        !is.null(internal_classes) &&
        !all(
          internal_classes %in%
          names(label_map)
        )
      ) {
        
        missing_labels <- internal_classes[
          !internal_classes %in%
            names(label_map)
        ]
        
        
        stop(
          paste0(
            "`label_map` does not contain mappings for: ",
            paste(
              missing_labels,
              collapse = ", "
            )
          )
        )
      }
      
      
      model_object$label_map <- label_map
      
      
      if (!is.null(internal_classes)) {
        
        model_object$classes <- unname(
          label_map[
            internal_classes
          ]
        )
      }
    }
    
    
    message(
      "Native AlgaeScan model loaded successfully."
    )
    
    
    if (!is.null(model_object$method)) {
      
      message(
        "Method: ",
        model_object$method
      )
    }
    
    
    return(
      model_object
    )
  }
  
  
  # ============================================================
  # 5. Check for legacy caret model
  # ============================================================
  
  if (
    !inherits(
      model_object,
      "train"
    )
  ) {
    
    stop(
      paste0(
        "The RDS file does not contain an `AlgaeScan_model` ",
        "or a supported legacy `caret::train` model."
      )
    )
  }
  
  
  message(
    "Legacy caret model detected."
  )
  
  
  # ============================================================
  # 6. Recover predictor names
  # ============================================================
  
  features <- NULL
  
  
  # ------------------------------------------------------------
  # Predictor names stored directly by caret
  # ------------------------------------------------------------
  
  if (!is.null(model_object$coefnames)) {
    
    features <- model_object$coefnames
  }
  
  
  # ------------------------------------------------------------
  # Some caret backends store predictor names inside finalModel.
  #
  # This fallback is important for legacy AlgaeScan Random Forest
  # models, but does not depend on the model method.
  # ------------------------------------------------------------
  
  if (
    is.null(features) &&
    !is.null(model_object$finalModel) &&
    !is.null(model_object$finalModel$xNames)
  ) {
    
    features <- model_object$finalModel$xNames
  }
  
  
  # ------------------------------------------------------------
  # Preprocessing statistics also contain predictor names
  # ------------------------------------------------------------
  
  if (
    is.null(features) &&
    !is.null(model_object$preProcess) &&
    !is.null(model_object$preProcess$mean)
  ) {
    
    features <- names(
      model_object$preProcess$mean
    )
  }
  
  
  # ------------------------------------------------------------
  # Final fallback: retained caret training data
  # ------------------------------------------------------------
  
  if (
    is.null(features) &&
    !is.null(model_object$trainingData)
  ) {
    
    features <- setdiff(
      colnames(
        model_object$trainingData
      ),
      ".outcome"
    )
  }
  
  
  if (!is.null(features)) {
    
    features <- as.character(
      features
    )
    
  } else {
    
    warning(
      paste0(
        "Predictor names could not be recovered from the legacy model. ",
        "`AlgaeScan_predict()` will use the available numeric columns ",
        "in the prediction data."
      ),
      call. = FALSE
    )
  }
  
  
  # ============================================================
  # 7. Recover model classes
  # ============================================================
  
  internal_classes <- NULL
  
  
  if (!is.null(model_object$levels)) {
    
    internal_classes <- as.character(
      model_object$levels
    )
  }
  
  
  if (
    is.null(internal_classes) &&
    !is.null(model_object$obsLevels)
  ) {
    
    internal_classes <- as.character(
      model_object$obsLevels
    )
  }
  
  
  if (
    is.null(internal_classes) &&
    !is.null(model_object$finalModel) &&
    !is.null(model_object$finalModel$classes)
  ) {
    
    internal_classes <- as.character(
      model_object$finalModel$classes
    )
  }
  
  
  # ============================================================
  # 8. Prepare class label map
  # ============================================================
  
  if (!is.null(label_map)) {
    
    if (
      !is.null(internal_classes) &&
      !all(
        internal_classes %in%
        names(label_map)
      )
    ) {
      
      missing_labels <- internal_classes[
        !internal_classes %in%
          names(label_map)
      ]
      
      
      stop(
        paste0(
          "`label_map` does not contain mappings for: ",
          paste(
            missing_labels,
            collapse = ", "
          )
        )
      )
    }
    
  } else if (!is.null(internal_classes)) {
    
    # Preserve legacy labels exactly when the user does not
    # request a different output representation.
    
    label_map <- stats::setNames(
      internal_classes,
      internal_classes
    )
  }
  
  
  # ============================================================
  # 9. Recover tuning parameters
  # ============================================================
  
  best_parameters <- list()
  
  
  if (!is.null(model_object$bestTune)) {
    
    best_parameters <- as.list(
      model_object$bestTune[
        1,
        ,
        drop = FALSE
      ]
    )
    
    
    # caret versions/models may use parameter names such as:
    #
    # .mtry
    # .k
    # .size
    #
    # Remove the leading dot generically instead of adding
    # model-specific handling.
    
    names(best_parameters) <- sub(
      "^\\.",
      "",
      names(best_parameters)
    )
  }
  
  
  # ============================================================
  # 10. Recover validation results
  # ============================================================
  
  validation_results <- NULL
  
  best_accuracy <- NA_real_
  
  
  if (!is.null(model_object$results)) {
    
    validation_results <- data.table::as.data.table(
      data.table::copy(
        model_object$results
      )
    )
    
    
    # Standardize tuning parameter names in the results table
    # using the same generic rule.
    
    parameter_names <- colnames(
      validation_results
    )
    
    
    cleaned_names <- sub(
      "^\\.",
      "",
      parameter_names
    )
    
    
    data.table::setnames(
      validation_results,
      parameter_names,
      cleaned_names
    )
    
    
    if (
      "Accuracy" %in%
      colnames(validation_results)
    ) {
      
      accuracy_values <- validation_results$Accuracy
      
      
      if (any(is.finite(accuracy_values))) {
        
        best_accuracy <- max(
          accuracy_values,
          na.rm = TRUE
        )
      }
    }
  }
  
  
  # ============================================================
  # 11. Determine output class labels
  # ============================================================
  
  if (
    !is.null(label_map) &&
    !is.null(internal_classes)
  ) {
    
    output_classes <- unname(
      label_map[
        internal_classes
      ]
    )
    
  } else {
    
    output_classes <- internal_classes
  }
  
  
  # ============================================================
  # 12. Build AlgaeScan wrapper
  # ============================================================
  
  output <- list(
    
    # Original fitted caret model
    model = model_object,
    
    # Algorithm used by caret
    method = model_object$method,
    
    # Target is optional because old caret objects may not
    # contain a meaningful original response-variable name
    target = target,
    
    # Predictor names recovered from the model
    features = features,
    
    # Class information
    classes = output_classes,
    label_map = label_map,
    
    # Generic caret tuning parameters
    best_parameters = best_parameters,
    
    # Available validation information
    best_accuracy = best_accuracy,
    validation_results = validation_results,
    
    # These fields were not part of the old AlgaeScan workflow
    best_repetition = NA_integer_,
    repetition_results = NULL,
    
    # Training-event metadata cannot reliably be reconstructed
    # from every legacy caret model
    training_events = NA_integer_,
    input_events = NA_integer_,
    input_class_counts = NULL,
    training_class_counts = NULL,
    
    # Mark this object as an imported legacy model
    settings = list(
      legacy_model = TRUE
    ),
    
    # Exact seed/software metadata are generally unavailable
    # for old caret objects
    reproducibility = list(
      source = "legacy caret model",
      software_versions = NULL
    ),
    
    training_time = NULL
  )
  
  
  class(
    output
  ) <- "AlgaeScan_model"
  
  
  # ============================================================
  # 13. Summary
  # ============================================================
  
  message(
    "Legacy caret model loaded successfully."
  )
  
  
  if (!is.null(model_object$method)) {
    
    message(
      "Method: ",
      model_object$method
    )
  }
  
  
  if (!is.null(features)) {
    
    message(
      "Predictor variables recovered: ",
      length(features)
    )
  }
  
  
  if (is.null(target)) {
    
    message(
      "Target variable: not specified"
    )
    
  } else {
    
    message(
      "Target variable: ",
      target
    )
  }
  
  
  return(
    output
  )
}


#' Save an AlgaeScan model
#'
#' Saves a trained AlgaeScan model to an RDS file.
#'
#' The complete `AlgaeScan_model` object is stored, including the
#' fitted model, predictor names, class labels, selected
#' hyperparameters, validation results, training settings and
#' reproducibility information.
#'
#' Models saved with this function can be restored using
#' [AlgaeScan_load_model()].
#'
#' @param model An object of class `AlgaeScan_model`, generated by
#'   [AlgaeScan_train()] or [AlgaeScan_load_model()].
#'
#' @param file Character string giving the path of the output RDS file.
#'
#' @param overwrite Logical. If `FALSE`, an existing file will not be
#'   overwritten. Default is `FALSE`.
#'
#' @param compress Compression passed to `saveRDS()`.
#'   Default is `TRUE`. Set to `FALSE` for faster saving and loading
#'   at the cost of a larger file.
#'
#' @return Invisibly returns the path of the saved model.
#'
#' @examples
#' \dontrun{
#'
#' AlgaeScan_save_model(
#'   model = model_A,
#'   file = "model_A.rds"
#' )
#'
#' }
#'
#' @export
AlgaeScan_save_model <- function(
    model,
    file,
    overwrite = FALSE,
    compress = TRUE
) {
  
  # ============================================================
  # 1. Check model
  # ============================================================
  
  if (!inherits(model, "AlgaeScan_model")) {
    
    stop(
      "`model` must be an object of class `AlgaeScan_model`."
    )
  }
  
  
  if (is.null(model$model)) {
    
    stop(
      "The AlgaeScan model does not contain a fitted model."
    )
  }
  
  
  # ============================================================
  # 2. Check output path
  # ============================================================
  
  if (
    !is.character(file) ||
    length(file) != 1 ||
    is.na(file) ||
    file == ""
  ) {
    
    stop(
      "`file` must be a single valid file path."
    )
  }
  
  
  file <- path.expand(
    file
  )
  
  
  output_dir <- dirname(
    file
  )
  
  
  if (!dir.exists(output_dir)) {
    
    stop(
      paste0(
        "Output directory does not exist:\n",
        output_dir
      )
    )
  }
  
  
  if (
    file.exists(file) &&
    !overwrite
  ) {
    
    stop(
      paste0(
        "File already exists:\n",
        file,
        "\nUse `overwrite = TRUE` to replace it."
      )
    )
  }
  
  
  # ============================================================
  # 3. Save first to a temporary file
  # ============================================================
  
  # Saving to a temporary file prevents an incomplete model file
  # from replacing an existing valid model if saveRDS() fails.
  
  temp_file <- tempfile(
    pattern = ".AlgaeScan_",
    tmpdir = output_dir,
    fileext = ".rds"
  )
  
  
  on.exit({
    
    if (file.exists(temp_file)) {
      
      unlink(
        temp_file
      )
    }
    
  }, add = TRUE)
  
  
  message(
    "Saving AlgaeScan model..."
  )
  
  
  saveRDS(
    object = model,
    file = temp_file,
    compress = compress
  )
  
  
  # ============================================================
  # 4. Move completed file to requested destination
  # ============================================================
  
  if (
    overwrite &&
    file.exists(file)
  ) {
    
    removed <- unlink(
      file
    )
    
    
    if (removed != 0) {
      
      stop(
        paste0(
          "Existing model file could not be replaced:\n",
          file
        )
      )
    }
  }
  
  
  moved <- file.rename(
    from = temp_file,
    to = file
  )
  
  
  if (!moved) {
    
    stop(
      paste0(
        "The model was created but could not be moved to:\n",
        file
      )
    )
  }
  
  
  # ============================================================
  # 5. Confirm saved model
  # ============================================================
  
  file_size_mb <- file.info(
    file
  )$size / 1024^2
  
  
  message(
    "Model saved successfully:"
  )
  
  message(
    file
  )
  
  message(
    "File size: ",
    round(
      file_size_mb,
      1
    ),
    " MB"
  )
  
  
  invisible(
    file
  )
}