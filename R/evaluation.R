#' Calculate F1 score
#'
#' Calculates binary or multiclass F1 scores for AlgaeScan predictions.
#'
#' For binary classification, the F1 score is calculated using
#' `MLmetrics::F1_Score()`.
#'
#' For multiclass classification, F1 is calculated independently for
#' each class from precision and recall obtained from
#' `caret::confusionMatrix()`. The macro-F1 is the mean of the
#' class-specific F1 scores.
#'
#' @param y_truth Vector containing reference class labels.
#' @param y_pred Vector containing predicted class labels.
#' @param pos_label Positive class for binary classification.
#'   Default is `"algae"`.
#' @param multiclass Logical. If `FALSE`, calculate binary F1.
#'   If `TRUE`, calculate class-specific and macro-F1 scores.
#'   Default is `FALSE`.
#'
#' @return For binary classification, a numeric F1 score.
#'   For multiclass classification, a list containing `macro_score`,
#'   `all_scores`, and `conf_matrix`.
#'
#' @export
AlgaeScan_f1 <- function(
    y_truth,
    y_pred,
    pos_label = "algae",
    multiclass = FALSE
) {
  
  # ------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------
  
  if (length(y_truth) != length(y_pred)) {
    stop("`y_truth` and `y_pred` must have the same length.")
  }
  
  if (length(y_truth) == 0) {
    stop("`y_truth` and `y_pred` contain no observations.")
  }
  
  
  y_truth <- as.character(y_truth)
  y_pred <- as.character(y_pred)
  
  
  # ============================================================
  # Binary F1
  # ============================================================
  
  if (!multiclass) {
    
    score <- MLmetrics::F1_Score(
      y_true = y_truth,
      y_pred = y_pred,
      positive = pos_label
    )
    
    score <- round(
      score,
      3
    )
    
    return(score)
  }
  
  
  # ============================================================
  # Multiclass F1
  # ============================================================
  
  y_truth <- as.factor(
    y_truth
  )
  
  y_pred <- as.factor(
    y_pred
  )
  
  
  # ------------------------------------------------------------
  # Align factor levels if required
  # ------------------------------------------------------------
  
  check_levels <- levels(y_truth) %in% levels(y_pred)
  
  
  if (any(!check_levels)) {
    
    all_levels <- union(
      levels(y_truth),
      levels(y_pred)
    )
    
    y_truth <- factor(
      y_truth,
      levels = all_levels
    )
    
    y_pred <- factor(
      y_pred,
      levels = all_levels
    )
  }
  
  
  # ------------------------------------------------------------
  # Confusion matrix
  # ------------------------------------------------------------
  
  conf_matrix <- caret::confusionMatrix(
    y_pred,
    y_truth
  )
  
  
  # ------------------------------------------------------------
  # F1 for each class
  # ------------------------------------------------------------
  
  precision <- conf_matrix$byClass[, "Precision"]
  
  recall <- conf_matrix$byClass[, "Recall"]
  
  
  f1_scores <- 2 * (precision * recall) /
    (precision + recall)
  
  
  macro_f1 <- mean(
    f1_scores,
    na.rm = TRUE
  )
  
  
  return(
    list(
      macro_score = macro_f1,
      all_scores = f1_scores,
      conf_matrix = conf_matrix
    )
  )
}


#' Plot a confusion matrix
#'
#' Creates a confusion matrix comparing reference and predicted
#' classification labels.
#'
#' Each cell shows the number of events and the percentage relative
#' to the corresponding reference class.
#'
#' @param y_truth Vector containing the reference class labels.
#'
#' @param y_pred Vector containing the predicted class labels.
#'
#' @param size_axis_text Numeric value controlling axis text size.
#'   Default is `18`.
#'
#' @param size_title_x Numeric value controlling the x-axis title size.
#'   Default is `20`.
#'
#' @param size_title_y Numeric value controlling the y-axis title size.
#'   Default is `20`.
#'
#' @param show_legend Logical. If `TRUE`, display the colour-scale legend.
#'   Default is `FALSE`.
#'
#' @return A `ggplot` object containing the confusion matrix.
#'
#' @examples
#' \dontrun{
#'
#' AlgaeScan_confusion_matrix(
#'   y_truth = prepared$test$type,
#'   y_pred = pred_A$prediction
#' )
#'
#' }
#'
#' @export
AlgaeScan_confusion_matrix <- function(
    y_truth,
    y_pred,
    size_axis_text = 18,
    size_title_x = 20,
    size_title_y = 20,
    show_legend = FALSE
) {
  
  # ------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------
  
  if (length(y_truth) != length(y_pred)) {
    stop(
      "`y_truth` and `y_pred` must have the same length."
    )
  }
  
  if (length(y_truth) == 0) {
    stop(
      "`y_truth` and `y_pred` contain no observations."
    )
  }
  
  if (any(is.na(y_truth)) || any(is.na(y_pred))) {
    stop(
      "`y_truth` and `y_pred` cannot contain missing values."
    )
  }
  
  
  y_truth <- as.character(
    y_truth
  )
  
  y_pred <- as.character(
    y_pred
  )
  
  
  # ============================================================
  # Build confusion matrix
  # ============================================================
  
  conf_matrix <- table(
    Predicted = y_pred,
    Reference = y_truth
  )
  
  
  conf_df <- as.data.frame(
    conf_matrix
  )
  
  
  # ------------------------------------------------------------
  # Calculate proportions within each reference class
  # ------------------------------------------------------------
  
  conf_df$Proportion <- mapply(
    function(freq, reference) {
      
      freq / sum(
        conf_df$Freq[
          conf_df$Reference == reference
        ]
      )
      
    },
    conf_df$Freq,
    conf_df$Reference
  )
  
  
  # ------------------------------------------------------------
  # Add count and percentage labels
  # ------------------------------------------------------------
  
  conf_df$Label <- sprintf(
    "%d\n(%.1f%%)",
    conf_df$Freq,
    conf_df$Proportion * 100
  )
  
  
  # ============================================================
  # Plot
  # ============================================================
  
  gg_plot <- ggplot2::ggplot(
    conf_df,
    ggplot2::aes(
      x = Reference,
      y = Predicted,
      fill = Proportion
    )
  ) +
    ggplot2::geom_tile(
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Label
      ),
      color = "white",
      size = 6
    ) +
    ggplot2::scale_fill_gradient(
      low = "#001F3F",
      high = "#E36414"
    ) +
    ggplot2::labs(
      title = NULL,
      x = "Reference",
      y = "Predicted"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        size = size_axis_text
      ),
      
      axis.title.x = ggplot2::element_text(
        size = size_title_x,
        face = "bold"
      ),
      
      axis.title.y = ggplot2::element_text(
        size = size_title_y,
        face = "bold"
      ),
      
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_blank(),
      
      panel.grid.major = ggplot2::element_blank(),
      
      legend.key.size = grid::unit(
        1,
        "cm"
      ),
      
      legend.title = ggplot2::element_text(
        size = 20
      ),
      
      legend.text = ggplot2::element_text(
        size = 20
      )
    )
  
  
  if (!show_legend) {
    
    gg_plot <- gg_plot +
      ggplot2::theme(
        legend.position = "none"
      )
  }
  
  
  return(gg_plot)
}