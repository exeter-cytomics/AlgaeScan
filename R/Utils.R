#' Collect repeated cross-validation results
#'
#' Combines repeated cross-validation results from multiple
#' `AlgaeScan_model` objects into a single long-format table.
#'
#' The function supports both supervised and unsupervised AlgaeScan
#' models. For supervised models, the cross-validation score is taken
#' from `best_accuracy`. For unsupervised models, it is taken from
#' `best_score`.
#'
#' Model names are taken from the names of the input list.
#'
#' @param models Named list of `AlgaeScan_model` objects containing
#'   `repetition_results`.
#'
#' @return A `data.table` containing:
#'
#' \itemize{
#'   \item `model`: model name taken from `names(models)`.
#'   \item `repetition`: training repetition number.
#'   \item `score`: best validation score for that repetition.
#'   \item `score_metric`: validation metric used by the model.
#' }
#'
#' @examples
#' \dontrun{
#'
#' cv_results <- AlgaeScan_collect_cv(
#'   list(
#'     RF = model_rf,
#'     KNN = model_knn,
#'     IF = model_if
#'   )
#' )
#'
#' }
#'
#' @export
AlgaeScan_collect_cv <- function(
    models
) {
  
  # ============================================================
  # 1. Check input
  # ============================================================
  
  if (!is.list(models)) {
    
    stop(
      "`models` must be a named list of AlgaeScan models."
    )
  }
  
  
  if (length(models) == 0) {
    
    stop(
      "`models` contains no models."
    )
  }
  
  
  if (
    is.null(names(models)) ||
    any(names(models) == "")
  ) {
    
    stop(
      "`models` must be a named list."
    )
  }
  
  
  # ============================================================
  # 2. Extract repetition results
  # ============================================================
  
  results <- lapply(
    names(models),
    function(model_name) {
      
      model <- models[[model_name]]
      
      
      if (is.null(model$repetition_results)) {
        
        stop(
          paste0(
            "Model `",
            model_name,
            "` does not contain `repetition_results`."
          )
        )
      }
      
      
      repetition_results <- data.table::as.data.table(
        data.table::copy(
          model$repetition_results
        )
      )
      
      
      if (!"repetition" %in% colnames(repetition_results)) {
        
        stop(
          paste0(
            "Model `",
            model_name,
            "` does not contain a `repetition` column."
          )
        )
      }
      
      
      # --------------------------------------------------------
      # Determine score column and metric
      # --------------------------------------------------------
      
      if ("best_accuracy" %in% colnames(repetition_results)) {
        
        score <- repetition_results$best_accuracy
        
        score_metric <- rep(
          "Accuracy",
          nrow(repetition_results)
        )
        
      } else if ("best_score" %in% colnames(repetition_results)) {
        
        score <- repetition_results$best_score
        
        
        if ("score_metric" %in% colnames(repetition_results)) {
          
          score_metric <- as.character(
            repetition_results$score_metric
          )
          
        } else if (!is.null(model$score_metric)) {
          
          score_metric <- rep(
            model$score_metric,
            nrow(repetition_results)
          )
          
        } else {
          
          score_metric <- rep(
            "Score",
            nrow(repetition_results)
          )
        }
        
      } else {
        
        stop(
          paste0(
            "Model `",
            model_name,
            "` contains neither `best_accuracy` nor `best_score`."
          )
        )
      }
      
      
      data.table::data.table(
        model = model_name,
        repetition = repetition_results$repetition,
        score = score,
        score_metric = score_metric
      )
    }
  )
  
  
  # ============================================================
  # 3. Combine models
  # ============================================================
  
  results <- data.table::rbindlist(
    results,
    use.names = TRUE
  )
  
  
  return(
    results
  )
}



#' Summarize repeated cross-validation scores
#'
#' Calculates summary statistics across cross-validation repetitions
#' for each model.
#'
#' The returned statistics include the mean, standard deviation,
#' minimum and maximum score across repetitions.
#'
#' @param data A `data.frame` or `data.table` containing at least
#'   `model` and `score` columns, such as the output returned by
#'   [AlgaeScan_collect_cv()].
#'
#' @param digits Number of decimal places used to round summary
#'   statistics. Default is `3`.
#'
#' @return A `data.table` containing one row per model with
#'   `mean_score`, `sd_score`, `min_score` and `max_score`.
#'
#' @examples
#' \dontrun{
#'
#' cv_summary <- AlgaeScan_summarize_cv(
#'   supp4_C_data
#' )
#'
#' }
#'
#' @export
AlgaeScan_summarize_cv <- function(
    data,
    digits = 3
) {
  
  # ============================================================
  # 1. Check input
  # ============================================================
  
  if (!is.data.frame(data)) {
    
    stop(
      "`data` must be a data.frame or data.table."
    )
  }
  
  
  required_cols <- c(
    "model",
    "score"
  )
  
  
  missing_cols <- required_cols[
    !required_cols %in% colnames(data)
  ]
  
  
  if (length(missing_cols) > 0) {
    
    stop(
      paste0(
        "Missing required column(s): ",
        paste(
          missing_cols,
          collapse = ", "
        )
      )
    )
  }
  
  
  if (
    !is.numeric(digits) ||
    length(digits) != 1 ||
    is.na(digits) ||
    digits < 0 ||
    digits %% 1 != 0
  ) {
    
    stop(
      "`digits` must be a single non-negative integer."
    )
  }
  
  
  data <- data.table::as.data.table(
    data.table::copy(data)
  )
  
  
  # ============================================================
  # 2. Summarize scores
  # ============================================================
  
  summary_data <- data[
    ,
    .(
      mean_score = mean(
        score,
        na.rm = TRUE
      ),
      sd_score = stats::sd(
        score,
        na.rm = TRUE
      ),
      min_score = min(
        score,
        na.rm = TRUE
      ),
      max_score = max(
        score,
        na.rm = TRUE
      )
    ),
    by = model
  ]
  
  
  # ============================================================
  # 3. Round output
  # ============================================================
  
  score_cols <- c(
    "mean_score",
    "sd_score",
    "min_score",
    "max_score"
  )
  
  
  summary_data[
    ,
    (score_cols) := lapply(
      .SD,
      round,
      digits = digits
    ),
    .SDcols = score_cols
  ]
  
  
  return(
    summary_data
  )
}