# ============================================================
# AlgaeScan model training
#
# Public function:
#   AlgaeScan_train()
#
# Internal model-specific functions:
#   .AlgaeScan_train_rf()
#   .AlgaeScan_train_knn()
#   .AlgaeScan_train_nnet()
#   .AlgaeScan_train_log_reg()
#   .AlgaeScan_train_nb()
#
# Internal shared training engine:
#   .AlgaeScan_train_caret()
#
# Internal reproducibility helpers:
#   .AlgaeScan_make_caret_seeds()
#   .AlgaeScan_package_version()
# ============================================================



#' Train an AlgaeScan model
#'
#' Trains a supervised, one-class or unsupervised machine-learning
#' model using spectral flow cytometry data.
#'
#' `AlgaeScan_train()` is the main user-facing training function. It
#' dispatches training to the appropriate internal model-specific
#' function according to the value supplied to `method`.
#'
#' Supported algorithms are:
#'
#' \itemize{
#'   \item `"rf"`: Random Forest
#'   \item `"knn"`: k-nearest neighbours
#'   \item `"nnet"`: neural network with one hidden layer
#'   \item `"log_reg"`: multinomial logistic regression
#'   \item `"nb"`: Naive Bayes
#'   \item `"svm"`: linear Support Vector Machine
#'   \item `"one_class_svm"`: One-Class Support Vector Machine
#'   \item `"gmm"`: Gaussian Mixture Model
#'   \item `"isoforest"`: Isolation Forest
#' }
#'
#' The first six methods are supervised classifiers.
#'
#' One-Class SVM and Isolation Forest are one-class anomaly-detection
#' methods. They are trained using events from a reference class,
#' specified using the model-specific `ref_label` argument, while
#' validation is performed using all classes.
#'
#' Gaussian Mixture Models are fitted without using the target labels.
#' Target labels are used only for cross-validation splitting and
#' evaluation. GMM performance is evaluated using the adjusted Rand
#' index rather than classification accuracy.
#'
#' Model-specific hyperparameters are supplied through `...`.
#'
#' For Random Forest:
#'
#' \code{
#' n_tree = c(5, 10, 15)
#' mtry = c(10, 20, 30, 40, 50, 60, 70)
#' }
#'
#' For KNN:
#'
#' \code{
#' k = c(3, 5, 7, 9, 11)
#' }
#'
#' For the neural network:
#'
#' \code{
#' hidden_units = c(3, 5, 7)
#' decay = c(0, 0.001, 0.01)
#' maxit = 100
#' maxnwts = 10000
#' }
#'
#' For multinomial logistic regression:
#'
#' \code{
#' decay = c(0.1, 0.01, 0.001)
#' maxnwts = 10000
#' }
#'
#' For Naive Bayes:
#'
#' \code{
#' laplace = c(0, 1)
#' usekernel = c(TRUE, FALSE)
#' adjust = c(1, 2)
#' }
#'
#' For the linear Support Vector Machine:
#'
#' \code{
#' cost = c(0.1, 1, 10)
#' }
#'
#' For One-Class SVM:
#'
#' \code{
#' ref_label = "algae"
#' nu = c(0.1, 0.2, 0.3, 0.4, 0.5)
#' }
#'
#' For Isolation Forest:
#'
#' \code{
#' ref_label = "algae"
#' contamination = c(0.1, 0.2, 0.3, 0.4, 0.5)
#' python_random_state = 42
#' }
#'
#' Isolation Forest is fitted using the Python scikit-learn
#' implementation distributed with AlgaeScan.
#'
#' The input data are copied internally and are never reordered by
#' class. When downsampling is requested, selected events retain
#' their original relative order.
#'
#' @param data A `data.frame` or `data.table` containing training data.
#'
#' @param target Character string giving the response or reference
#'   variable. For supervised models this is the classification
#'   target. For unsupervised models it is used for validation.
#'
#' @param method Training algorithm. One of `"rf"`, `"knn"`, `"nnet"`,
#'   `"log_reg"`, `"nb"`, `"svm"`, `"one_class_svm"`, `"gmm"` or
#'   `"isoforest"`. Default is `"rf"`.
#'
#' @param feature_cols Optional character vector containing predictor
#'   columns. If `NULL`, all numeric columns except `target` are used.
#'
#' @param downsample Optional integer giving the maximum number of
#'   events retained per group. If `NULL`, all events are used.
#'
#' @param downsample_by Character vector containing the column or
#'   columns defining downsampling groups. Default is `"Species"`.
#'
#' @param cv_folds Number of cross-validation folds. Default is `5`.
#'
#' @param n_repeats Number of complete training repetitions.
#'   Default is `1`.
#'
#' @param method_control Validation method. `"cv"` performs
#'   cross-validation. `"oob"` performs out-of-bag validation and is
#'   currently available only for Random Forest.
#'
#' @param preprocess Optional preprocessing applied during training.
#'   Default is `c("center", "scale")`. For supervised `caret` models,
#'   preprocessing is handled by `caret::train()`. For unsupervised
#'   models, the fitted preprocessing object is stored with the
#'   AlgaeScan model and reused during prediction.
#'
#' @param n_cores Number of CPU cores requested during training.
#'   Default is `1`. The current unsupervised training engine runs
#'   cross-validation folds sequentially, so this argument does not
#'   currently parallelize One-Class SVM, GMM or Isolation Forest.
#'
#' @param seed Integer controlling reproducible training.
#'   Default is `123`.
#'
#' @param ... Model-specific training arguments.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @examples
#' \dontrun{
#'
#' # Random Forest
#' model_rf <- AlgaeScan_train(
#'   data = train_data,
#'   target = "type",
#'   method = "rf",
#'   n_tree = c(5, 10, 15),
#'   mtry = c(10, 20, 30, 40, 50, 60, 70),
#'   cv_folds = 5,
#'   n_cores = 4,
#'   seed = 123
#' )
#'
#' # One-Class SVM
#' model_ocsvm <- AlgaeScan_train(
#'   data = train_data,
#'   target = "type",
#'   method = "one_class_svm",
#'   ref_label = "algae",
#'   nu = c(0.1, 0.2, 0.3, 0.4, 0.5),
#'   seed = 123
#' )
#'
#' # Gaussian Mixture Model
#' model_gmm <- AlgaeScan_train(
#'   data = train_data,
#'   target = "type",
#'   method = "gmm",
#'   preprocess = NULL,
#'   seed = 123
#' )
#'
#' # Isolation Forest
#' model_if <- AlgaeScan_train(
#'   data = train_data,
#'   target = "type",
#'   method = "isoforest",
#'   ref_label = "algae",
#'   contamination = c(0.1, 0.2, 0.3, 0.4, 0.5),
#'   seed = 123
#' )
#'
#' }
#'
#' @export
AlgaeScan_train <- function(
    data,
    target,
    method = "rf",
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    ...
) {
  
  # ------------------------------------------------------------
  # Select training algorithm
  # ------------------------------------------------------------
  
  method <- match.arg(
    method,
    c(
      "rf",
      "knn",
      "nnet",
      "log_reg",
      "nb",
      "svm",
      "one_class_svm",
      "gmm",
      "isoforest"
    )
  )
  
  # ------------------------------------------------------------
  # Arguments shared by all training methods
  # ------------------------------------------------------------
  
  common_args <- list(
    data = data,
    target = target,
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed
  )
  
  
  # Model-specific arguments supplied through ...
  model_args <- list(...)
  
  
  # ------------------------------------------------------------
  # Dispatch to model-specific training function
  # ------------------------------------------------------------
  
  training_function <- switch(
    method,
    
    rf = .AlgaeScan_train_rf,
    
    knn = .AlgaeScan_train_knn,
    
    nnet = .AlgaeScan_train_nnet,
    
    log_reg = .AlgaeScan_train_log_reg,
    
    nb = .AlgaeScan_train_nb,
    
    svm = .AlgaeScan_train_svm,
    
    one_class_svm = .AlgaeScan_train_one_class_svm,
    
    gmm = .AlgaeScan_train_gmm,
    
    isoforest = .AlgaeScan_train_isoforest
  )
  
  
  do.call(
    training_function,
    c(
      common_args,
      model_args
    )
  )
}


#' Train an AlgaeScan Random Forest model
#'
#' Internal Random Forest training function used by
#' [AlgaeScan_train()].
#'
#' Random Forest models are tuned over combinations of `mtry` and
#' `n_tree`. Each value of `n_tree` is treated as a separate model
#' configuration, while `mtry` is tuned using `caret::train()`.
#'
#' @param data Training data.
#' @param target Response variable.
#' @param feature_cols Predictor columns.
#' @param downsample Optional number of events retained per group.
#' @param downsample_by Downsampling grouping variables.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method.
#' @param preprocess Preprocessing passed to `caret::train()`.
#' @param n_cores Number of CPU cores.
#' @param seed Reproducibility seed.
#' @param n_tree Integer vector containing numbers of trees.
#' @param mtry Integer vector containing candidate mtry values.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_rf <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    n_tree = c(5, 10, 15),
    mtry = c(10, 20, 30, 40, 50, 60, 70)
) {
  
  # ============================================================
  # 1. Check Random Forest parameters
  # ============================================================
  
  if (
    !is.numeric(n_tree) ||
    length(n_tree) == 0 ||
    any(n_tree <= 0) ||
    any(n_tree %% 1 != 0)
  ) {
    stop(
      "`n_tree` must contain positive integers."
    )
  }
  
  
  if (
    !is.numeric(mtry) ||
    length(mtry) == 0 ||
    any(mtry <= 0) ||
    any(mtry %% 1 != 0)
  ) {
    stop(
      "`mtry` must contain positive integers."
    )
  }
  
  
  n_tree <- unique(
    as.integer(n_tree)
  )
  
  
  mtry <- unique(
    as.integer(mtry)
  )
  
  
  # ============================================================
  # 2. Build Random Forest configurations
  # ============================================================
  
  # caret tunes mtry directly.
  #
  # The number of trees is passed directly to randomForest and
  # therefore requires one caret training run for each n_tree value.
  
  configurations <- lapply(
    n_tree,
    function(tree_value) {
      
      list(
        
        label = paste0(
          "n_tree = ",
          tree_value
        ),
        
        tune_grid = data.frame(
          mtry = mtry
        ),
        
        train_args = list(
          ntree = tree_value
        ),
        
        fixed_parameters = list(
          n_tree = tree_value
        )
      )
    }
  )
  
  
  # ============================================================
  # 3. Train model
  # ============================================================
  
  output <- .AlgaeScan_train_caret(
    data = data,
    target = target,
    method = "rf",
    method_name = "Random Forest",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    configurations = configurations,
    tuning_settings = list(
      n_tree = n_tree,
      mtry = mtry
    ),
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    backend_package = "randomForest"
  )
  
  
  # ============================================================
  # 4. Add convenient Random Forest fields
  # ============================================================
  
  output$best_n_tree <- output$best_parameters$n_tree
  
  output$best_mtry <- output$best_parameters$mtry
  
  
  return(
    output
  )
}



#' Train an AlgaeScan KNN model
#'
#' Internal k-nearest neighbours training function used by
#' [AlgaeScan_train()].
#'
#' Candidate values of `k` are evaluated using cross-validation.
#'
#' @param data Training data.
#' @param target Response variable.
#' @param feature_cols Predictor columns.
#' @param downsample Optional number of events retained per group.
#' @param downsample_by Downsampling grouping variables.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method.
#' @param preprocess Preprocessing passed to `caret::train()`.
#' @param n_cores Number of CPU cores.
#' @param seed Reproducibility seed.
#' @param k Integer vector containing candidate numbers of neighbours.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_knn <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    k = c(3, 5, 7, 9, 11)
) {
  
  # ============================================================
  # 1. Check KNN parameters
  # ============================================================
  
  if (method_control == "oob") {
    
    stop(
      "`method_control = \"oob\"` is not available for KNN."
    )
  }
  
  
  if (
    !is.numeric(k) ||
    length(k) == 0 ||
    any(k <= 0) ||
    any(k %% 1 != 0)
  ) {
    
    stop(
      "`k` must contain positive integers."
    )
  }
  
  
  k <- unique(
    as.integer(k)
  )
  
  
  # ============================================================
  # 2. Build KNN configuration
  # ============================================================
  
  configurations <- list(
    
    list(
      
      label = "KNN tuning",
      
      tune_grid = data.frame(
        k = k
      ),
      
      train_args = list(),
      
      fixed_parameters = list()
    )
  )
  
  
  # ============================================================
  # 3. Train model
  # ============================================================
  
  output <- .AlgaeScan_train_caret(
    data = data,
    target = target,
    method = "knn",
    method_name = "k-nearest neighbours",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    configurations = configurations,
    tuning_settings = list(
      k = k
    ),
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed
  )
  
  
  # ============================================================
  # 4. Add convenient KNN field
  # ============================================================
  
  output$best_k <- output$best_parameters$k
  
  
  return(
    output
  )
}



#' Train an AlgaeScan neural-network model
#'
#' Internal neural-network training function used by
#' [AlgaeScan_train()].
#'
#' The model uses `caret` method `"nnet"` and therefore contains a
#' single hidden layer. `hidden_units` controls the number of neurons
#' in that hidden layer, while `decay` controls weight decay.
#'
#' @param data Training data.
#' @param target Response variable.
#' @param feature_cols Predictor columns.
#' @param downsample Optional number of events retained per group.
#' @param downsample_by Downsampling grouping variables.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method.
#' @param preprocess Preprocessing passed to `caret::train()`.
#' @param n_cores Number of CPU cores.
#' @param seed Reproducibility seed.
#' @param hidden_units Integer vector giving numbers of neurons in the
#'   hidden layer.
#' @param decay Numeric vector giving candidate weight-decay values.
#' @param maxit Maximum number of neural-network iterations.
#' @param maxnwts Maximum number of neural-network weights.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_nnet <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    hidden_units = c(3, 5, 7),
    decay = c(0, 0.001, 0.01),
    maxit = 100,
    maxnwts = 10000
) {
  
  # ============================================================
  # 1. Check neural-network parameters
  # ============================================================
  
  if (method_control == "oob") {
    
    stop(
      "`method_control = \"oob\"` is not available for neural networks."
    )
  }
  
  
  if (
    !is.numeric(hidden_units) ||
    length(hidden_units) == 0 ||
    any(hidden_units <= 0) ||
    any(hidden_units %% 1 != 0)
  ) {
    
    stop(
      "`hidden_units` must contain positive integers."
    )
  }
  
  
  if (
    !is.numeric(decay) ||
    length(decay) == 0 ||
    any(decay < 0)
  ) {
    
    stop(
      "`decay` must contain non-negative numeric values."
    )
  }
  
  
  if (
    !is.numeric(maxit) ||
    length(maxit) != 1 ||
    maxit <= 0 ||
    maxit %% 1 != 0
  ) {
    
    stop(
      "`maxit` must be a single positive integer."
    )
  }
  
  
  if (
    !is.numeric(maxnwts) ||
    length(maxnwts) != 1 ||
    maxnwts <= 0 ||
    maxnwts %% 1 != 0
  ) {
    
    stop(
      "`maxnwts` must be a single positive integer."
    )
  }
  
  
  hidden_units <- unique(
    as.integer(hidden_units)
  )
  
  
  decay <- unique(
    decay
  )
  
  
  maxit <- as.integer(
    maxit
  )
  
  
  maxnwts <- as.integer(
    maxnwts
  )
  
  
  # ============================================================
  # 2. Build neural-network configuration
  # ============================================================
  
  # caret uses the name "size" for the number of hidden-layer
  # neurons. AlgaeScan exposes this more clearly as "hidden_units".
  
  configurations <- list(
    
    list(
      
      label = "Neural-network tuning",
      
      tune_grid = expand.grid(
        size = hidden_units,
        decay = decay
      ),
      
      train_args = list(
        trace = FALSE,
        maxit = maxit,
        MaxNWts = maxnwts
      ),
      
      fixed_parameters = list()
    )
  )
  
  
  # ============================================================
  # 3. Train model
  # ============================================================
  
  output <- .AlgaeScan_train_caret(
    data = data,
    target = target,
    method = "nnet",
    method_name = "Neural network",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    configurations = configurations,
    tuning_settings = list(
      hidden_units = hidden_units,
      decay = decay,
      maxit = maxit,
      maxnwts = maxnwts
    ),
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    parameter_map = c(
      "size" = "hidden_units"
    ),
    backend_package = "nnet"
  )
  
  
  # ============================================================
  # 4. Add convenient neural-network fields
  # ============================================================
  
  output$best_hidden_units <-
    output$best_parameters$hidden_units
  
  output$best_decay <-
    output$best_parameters$decay
  
  
  return(
    output
  )
}

#' Train an AlgaeScan multinomial logistic regression model
#'
#' Internal multinomial logistic regression training function used by
#' [AlgaeScan_train()].
#'
#' The model uses the `caret` method `"multinom"`. Candidate values
#' of `decay` are evaluated using cross-validation.
#'
#' @param data Training data.
#' @param target Response variable.
#' @param feature_cols Predictor columns.
#' @param downsample Optional number of events retained per group.
#' @param downsample_by Downsampling grouping variables.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method.
#' @param preprocess Preprocessing passed to `caret::train()`.
#' @param n_cores Number of CPU cores.
#' @param seed Reproducibility seed.
#' @param decay Numeric vector containing candidate weight-decay values.
#' @param maxnwts Maximum number of model weights.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_log_reg <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    decay = c(0.1, 0.01, 0.001),
    maxnwts = 10000
) {
  
  # ============================================================
  # 1. Check multinomial logistic regression parameters
  # ============================================================
  
  if (method_control == "oob") {
    
    stop(
      paste0(
        "`method_control = \"oob\"` is not available for ",
        "multinomial logistic regression."
      )
    )
  }
  
  
  if (
    !is.numeric(decay) ||
    length(decay) == 0 ||
    any(is.na(decay)) ||
    any(decay < 0)
  ) {
    
    stop(
      "`decay` must contain non-negative numeric values."
    )
  }
  
  
  if (
    !is.numeric(maxnwts) ||
    length(maxnwts) != 1 ||
    is.na(maxnwts) ||
    maxnwts <= 0 ||
    maxnwts %% 1 != 0
  ) {
    
    stop(
      "`maxnwts` must be a single positive integer."
    )
  }
  
  
  decay <- unique(
    decay
  )
  
  
  maxnwts <- as.integer(
    maxnwts
  )
  
  
  # ============================================================
  # 2. Build multinomial logistic regression configuration
  # ============================================================
  
  configurations <- list(
    
    list(
      
      label = "Multinomial logistic regression tuning",
      
      tune_grid = data.frame(
        decay = decay
      ),
      
      train_args = list(
        trace = FALSE,
        MaxNWts = maxnwts
      ),
      
      fixed_parameters = list()
    )
  )
  
  
  # ============================================================
  # 3. Train model
  # ============================================================
  
  output <- .AlgaeScan_train_caret(
    data = data,
    target = target,
    method = "multinom",
    method_name = "Multinomial logistic regression",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    configurations = configurations,
    tuning_settings = list(
      decay = decay,
      maxnwts = maxnwts
    ),
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    backend_package = "nnet"
  )
  
  
  # ============================================================
  # 4. Add convenient multinomial logistic regression field
  # ============================================================
  
  output$best_decay <-
    output$best_parameters$decay
  
  
  return(
    output
  )
}

#' Train an AlgaeScan Naive Bayes model
#'
#' Internal Naive Bayes training function used by
#' [AlgaeScan_train()].
#'
#' The model uses the `caret` method `"naive_bayes"`. Candidate
#' combinations of `laplace`, `usekernel` and `adjust` are evaluated
#' using cross-validation.
#'
#' @param data Training data.
#' @param target Response variable.
#' @param feature_cols Predictor columns.
#' @param downsample Optional number of events retained per group.
#' @param downsample_by Downsampling grouping variables.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method.
#' @param preprocess Preprocessing passed to `caret::train()`.
#' @param n_cores Number of CPU cores.
#' @param seed Reproducibility seed.
#' @param laplace Numeric vector containing candidate Laplace
#'   correction values.
#' @param usekernel Logical vector indicating whether kernel density
#'   estimation should be used.
#' @param adjust Numeric vector containing candidate bandwidth
#'   adjustment values.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_nb <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    laplace = c(0, 1),
    usekernel = c(TRUE, FALSE),
    adjust = c(1, 2)
) {
  
  # ============================================================
  # 1. Check Naive Bayes parameters
  # ============================================================
  
  if (method_control == "oob") {
    
    stop(
      "`method_control = \"oob\"` is not available for Naive Bayes."
    )
  }
  
  
  if (
    !is.numeric(laplace) ||
    length(laplace) == 0 ||
    any(is.na(laplace)) ||
    any(laplace < 0)
  ) {
    
    stop(
      "`laplace` must contain non-negative numeric values."
    )
  }
  
  
  if (
    !is.logical(usekernel) ||
    length(usekernel) == 0 ||
    any(is.na(usekernel))
  ) {
    
    stop(
      "`usekernel` must contain logical values."
    )
  }
  
  
  if (
    !is.numeric(adjust) ||
    length(adjust) == 0 ||
    any(is.na(adjust)) ||
    any(adjust <= 0)
  ) {
    
    stop(
      "`adjust` must contain positive numeric values."
    )
  }
  
  
  laplace <- unique(
    laplace
  )
  
  
  usekernel <- unique(
    usekernel
  )
  
  
  adjust <- unique(
    adjust
  )
  
  
  # ============================================================
  # 2. Build Naive Bayes configuration
  # ============================================================
  
  configurations <- list(
    
    list(
      
      label = "Naive Bayes tuning",
      
      tune_grid = expand.grid(
        laplace = laplace,
        usekernel = usekernel,
        adjust = adjust
      ),
      
      train_args = list(),
      
      fixed_parameters = list()
    )
  )
  
  
  # ============================================================
  # 3. Train model
  # ============================================================
  
  output <- .AlgaeScan_train_caret(
    data = data,
    target = target,
    method = "naive_bayes",
    method_name = "Naive Bayes",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    configurations = configurations,
    tuning_settings = list(
      laplace = laplace,
      usekernel = usekernel,
      adjust = adjust
    ),
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    backend_package = "naivebayes"
  )
  
  
  # ============================================================
  # 4. Add convenient Naive Bayes fields
  # ============================================================
  
  output$best_laplace <-
    output$best_parameters$laplace
  
  output$best_usekernel <-
    output$best_parameters$usekernel
  
  output$best_adjust <-
    output$best_parameters$adjust
  
  
  return(
    output
  )
}

#' Train an AlgaeScan linear SVM model
#'
#' Internal linear Support Vector Machine training function used by
#' [AlgaeScan_train()].
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_svm <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    cost = c(0.1, 1, 10)
) {
  
  # ============================================================
  # 1. Check SVM parameters
  # ============================================================
  
  if (method_control == "oob") {
    
    stop(
      "`method_control = \"oob\"` is not available for SVM."
    )
  }
  
  
  if (
    !is.numeric(cost) ||
    length(cost) == 0 ||
    any(is.na(cost)) ||
    any(cost <= 0)
  ) {
    
    stop(
      "`cost` must contain positive numeric values."
    )
  }
  
  
  cost <- unique(
    cost
  )
  
  
  # ============================================================
  # 2. Build SVM configuration
  # ============================================================
  
  configurations <- list(
    
    list(
      
      label = "Linear SVM tuning",
      
      tune_grid = data.frame(
        C = cost
      ),
      
      train_args = list(),
      
      fixed_parameters = list()
    )
  )
  
  
  # ============================================================
  # 3. Train model
  # ============================================================
  
  output <- .AlgaeScan_train_caret(
    data = data,
    target = target,
    method = "svmLinear",
    method_name = "Support Vector Machine",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    configurations = configurations,
    tuning_settings = list(
      cost = cost
    ),
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    method_control = method_control,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    parameter_map = c(
      "C" = "cost"
    ),
    backend_package = "kernlab"
  )
  
  
  output$best_cost <-
    output$best_parameters$cost
  
  
  return(
    output
  )
}

#' Train an AlgaeScan One-Class Support Vector Machine model
#'
#' Internal One-Class Support Vector Machine training function used by
#' [AlgaeScan_train()].
#'
#' The model learns the distribution of a single reference class,
#' defined by `ref_label`. For the AlgaeScan Model A workflow, the
#' reference class is `"algae"` and all other events are treated as
#' outliers.
#'
#' The model is fitted using `e1071::svm()` with
#' `type = "one-classification"` and a radial kernel.
#'
#' For reproduction of the original Supplementary Figure 4 workflow,
#' candidate `nu` values are assigned sequentially to the
#' cross-validation folds.
#'
#' @param data Training data as a `data.frame` or `data.table`.
#' @param target Character string giving the reference-label column.
#' @param feature_cols Optional character vector containing predictor
#'   columns. If `NULL`, numeric predictor columns are selected
#'   automatically.
#' @param downsample Optional maximum number of events retained per
#'   downsampling group.
#' @param downsample_by Character vector defining the downsampling
#'   groups.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method. Only `"cv"` is supported.
#' @param preprocess Preprocessing applied before cross-validation.
#'   Default is `c("center", "scale")`.
#' @param n_cores Number of CPU cores requested by the common training
#'   interface.
#' @param seed Integer controlling reproducible downsampling and
#'   cross-validation.
#' @param ref_label Character string defining the reference or normal
#'   class. Default is `"algae"`.
#' @param nu Numeric vector containing one-class SVM `nu` values.
#'   For reproduction of the original analysis, use
#'   `c(0.1, 0.2, 0.3, 0.4, 0.5)`.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_one_class_svm <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    ref_label = "algae",
    nu = c(0.1, 0.2, 0.3, 0.4, 0.5)
) {
  
  if (method_control != "cv") {
    
    stop(
      "One-Class SVM currently supports only cross-validation."
    )
  }
  
  
  if (
    !is.numeric(nu) ||
    length(nu) == 0 ||
    any(is.na(nu)) ||
    any(nu <= 0) ||
    any(nu > 1)
  ) {
    
    stop(
      "`nu` must contain values greater than 0 and at most 1."
    )
  }
  
  
  output <- .AlgaeScan_train_unsupervised(
    data = data,
    target = target,
    method = "one_class_svm",
    method_name = "One-Class SVM",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    ref_label = ref_label,
    fold_parameter_name = "nu",
    fold_parameter_values = nu,
    score_metric = "Accuracy"
  )
  
  
  output$best_nu <-
    output$best_parameters$nu
  
  output$best_accuracy <-
    output$best_score
  
  
  return(
    output
  )
}


#' Train an AlgaeScan Gaussian Mixture Model
#'
#' Internal Gaussian Mixture Model training function used by
#' [AlgaeScan_train()].
#'
#' Gaussian mixture models are fitted using `mclust::Mclust()`.
#' The number of mixture components and covariance structure are
#' selected by `mclust` using the Bayesian information criterion.
#'
#' The model itself is fitted without class labels. Reference labels
#' supplied through `target` are used only to construct
#' cross-validation folds and to evaluate the resulting clustering
#' using the adjusted Rand index.
#'
#' To reproduce the original AlgaeScan Supplementary Figure 4
#' analysis, use `preprocess = NULL`.
#'
#' @param data Training data as a `data.frame` or `data.table`.
#' @param target Character string containing reference labels used
#'   for validation only.
#' @param feature_cols Optional character vector containing predictor
#'   columns.
#' @param downsample Optional maximum number of events retained per
#'   downsampling group.
#' @param downsample_by Character vector defining the downsampling
#'   groups.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method. Only `"cv"` is supported.
#' @param preprocess Optional preprocessing. Use `NULL` to reproduce
#'   the original analysis.
#' @param n_cores Number of CPU cores requested by the common training
#'   interface.
#' @param seed Integer controlling reproducible downsampling and
#'   cross-validation.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_gmm <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = NULL,
    n_cores = 1,
    seed = 123
) {
  
  if (method_control != "cv") {
    
    stop(
      "GMM currently supports only cross-validation."
    )
  }
  
  
  output <- .AlgaeScan_train_unsupervised(
    data = data,
    target = target,
    method = "gmm",
    method_name = "Gaussian Mixture Model",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    ref_label = NULL,
    fold_parameter_name = NULL,
    fold_parameter_values = NULL,
    score_metric = "Adjusted Rand Index"
  )
  
  
  output$best_ari <-
    output$best_score
  
  
  return(
    output
  )
}

#' Train an AlgaeScan Isolation Forest model
#'
#' Internal Isolation Forest training function used by
#' [AlgaeScan_train()].
#'
#' Isolation Forest is used as a one-class anomaly-detection method.
#' The model is trained using only events belonging to the reference
#' class, while events from all classes are retained in the validation
#' folds.
#'
#' The Isolation Forest model itself is fitted in Python using the
#' scikit-learn implementation distributed with AlgaeScan. R handles
#' downsampling, preprocessing, cross-validation, repeated training
#' and model selection.
#'
#' For reproduction of the original AlgaeScan model-comparison
#' simulations, the contamination values
#' `c(0.1, 0.2, 0.3, 0.4, 0.5)` are assigned sequentially to the
#' five cross-validation folds.
#'
#' @param data Training data as a `data.frame` or `data.table`.
#'
#' @param target Character string giving the column containing the
#'   reference labels.
#'
#' @param feature_cols Optional character vector containing predictor
#'   columns. If `NULL`, numeric predictor columns are selected
#'   automatically.
#'
#' @param downsample Optional integer giving the maximum number of
#'   events retained per downsampling group. If `NULL`, all events
#'   are used.
#'
#' @param downsample_by Character vector containing the column or
#'   columns defining downsampling groups.
#'
#' @param cv_folds Number of cross-validation folds.
#'   Default is `5`.
#'
#' @param n_repeats Number of complete training repetitions.
#'   Default is `1`.
#'
#' @param method_control Validation method. Only `"cv"` is currently
#'   supported for Isolation Forest.
#'
#' @param preprocess Optional preprocessing applied before model
#'   training. Default is `c("center", "scale")`.
#'
#' @param n_cores Number of CPU cores requested by the common
#'   AlgaeScan training interface.
#'
#' @param seed Integer controlling reproducible downsampling and
#'   cross-validation. Default is `123`.
#'
#' @param ref_label Character string defining the normal or reference
#'   class used to train the Isolation Forest.
#'   Default is `"algae"`.
#'
#' @param contamination Numeric vector containing Isolation Forest
#'   contamination values. For reproduction of the original analysis,
#'   use `c(0.1, 0.2, 0.3, 0.4, 0.5)`.
#'
#' @param python_random_state Integer passed to scikit-learn's
#'   Isolation Forest as `random_state`.
#'   Default is `42`.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_isoforest <- function(
    data,
    target,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    ref_label = "algae",
    contamination = c(0.1, 0.2, 0.3, 0.4, 0.5),
    python_random_state = 42
) {
  
  # ============================================================
  # 1. Check Isolation Forest parameters
  # ============================================================
  
  if (method_control != "cv") {
    
    stop(
      "`method_control = \"oob\"` is not available for Isolation Forest."
    )
  }
  
  
  if (
    !is.numeric(contamination) ||
    length(contamination) == 0 ||
    any(is.na(contamination)) ||
    any(contamination <= 0) ||
    any(contamination > 0.5)
  ) {
    
    stop(
      "`contamination` must contain values greater than 0 and at most 0.5."
    )
  }
  
  
  if (
    length(contamination) != 1 &&
    length(contamination) != cv_folds
  ) {
    
    stop(
      paste0(
        "`contamination` must contain either one value or ",
        "one value per cross-validation fold."
      )
    )
  }
  
  
  if (
    !is.character(ref_label) ||
    length(ref_label) != 1 ||
    is.na(ref_label) ||
    ref_label == ""
  ) {
    
    stop(
      "`ref_label` must be a single class label."
    )
  }
  
  
  if (
    !is.numeric(python_random_state) ||
    length(python_random_state) != 1 ||
    is.na(python_random_state)
  ) {
    
    stop(
      "`python_random_state` must be a single integer."
    )
  }
  
  
  python_random_state <- as.integer(
    python_random_state
  )
  
  
  # ============================================================
  # 2. Train Isolation Forest
  # ============================================================
  
  output <- .AlgaeScan_train_unsupervised(
    data = data,
    target = target,
    method = "isoforest",
    method_name = "Isolation Forest",
    feature_cols = feature_cols,
    downsample = downsample,
    downsample_by = downsample_by,
    cv_folds = cv_folds,
    n_repeats = n_repeats,
    preprocess = preprocess,
    n_cores = n_cores,
    seed = seed,
    ref_label = ref_label,
    fold_parameter_name = "contamination",
    fold_parameter_values = contamination,
    score_metric = "Accuracy",
    python_random_state = python_random_state
  )
  
  
  # ============================================================
  # 3. Add convenient Isolation Forest fields
  # ============================================================
  
  output$best_contamination <-
    output$best_parameters$contamination
  
  output$best_accuracy <-
    output$best_score
  
  
  return(
    output
  )
}


#' Load the AlgaeScan Python Isolation Forest implementation
#'
#' Internal helper that loads the Python functions distributed with
#' AlgaeScan for fitting and predicting with scikit-learn Isolation
#' Forest models.
#'
#' @return An environment containing the Python functions.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_load_isoforest_python <- function() {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    
    stop(
      "Package `reticulate` is required for Isolation Forest."
    )
  }
  
  
  if (!reticulate::py_module_available("sklearn")) {
    
    stop(
      paste0(
        "Python module `sklearn` was not found. ",
        "Install scikit-learn in the Python environment used by reticulate."
      )
    )
  }
  
  
  python_file <- system.file(
    "python",
    "isolationforest_code.py",
    package = "AlgaeScan"
  )
  
  
  if (python_file == "") {
    
    stop(
      "AlgaeScan Isolation Forest Python code was not found."
    )
  }
  
  
  python_environment <- new.env(
    parent = baseenv()
  )
  
  
  reticulate::source_python(
    python_file,
    envir = python_environment
  )
  
  
  return(
    python_environment
  )
}

#' Shared caret training engine for AlgaeScan
#'
#' Internal training engine shared by the model-specific AlgaeScan
#' training functions.
#'
#' This function handles functionality that should behave identically
#' across supervised algorithms, including:
#'
#' \itemize{
#'   \item predictor selection
#'   \item optional downsampling
#'   \item preservation of input row order
#'   \item repeated training
#'   \item explicit cross-validation folds
#'   \item deterministic random seeds
#'   \item parallel processing
#'   \item model selection
#'   \item reproducibility metadata
#' }
#'
#' Model-specific tuning grids and additional arguments are supplied
#' through `configurations`.
#'
#' @param data Training data.
#' @param target Response variable.
#' @param method caret model identifier for the training method
#' @param method_name Human-readable model name for the training method
#' @param feature_cols Predictor columns.
#' @param downsample Optional number of events retained per group.
#' @param downsample_by Downsampling grouping variables.
#' @param configurations List of model configurations.
#' @param tuning_settings Model-specific settings stored in the output.
#' @param cv_folds Number of cross-validation folds.
#' @param n_repeats Number of complete training repetitions.
#' @param method_control Validation method.
#' @param preprocess Preprocessing passed to `caret::train()`.
#' @param n_cores Number of CPU cores.
#' @param seed Reproducibility seed.
#' @param parameter_map Optional named character vector used to rename
#'   caret parameter names in the AlgaeScan output.
#' @param backend_package Optional package name used by the model.
#'
#' @return An object of class `AlgaeScan_model`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_caret <- function(
    data,
    target,
    method,
    method_name,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    configurations,
    tuning_settings,
    cv_folds = 5,
    n_repeats = 1,
    method_control = "cv",
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    parameter_map = NULL,
    backend_package = NULL
) {
  
  # ============================================================
  # 1. General input checks
  # ============================================================
  
  if (!is.data.frame(data)) {
    
    stop(
      "`data` must be a data.frame or data.table."
    )
  }
  
  
  if (nrow(data) == 0) {
    
    stop(
      "`data` contains no events."
    )
  }
  
  
  if (
    !is.character(target) ||
    length(target) != 1 ||
    is.na(target) ||
    target == ""
  ) {
    
    stop(
      "`target` must be a single column name."
    )
  }
  
  
  if (!target %in% colnames(data)) {
    
    stop(
      paste0(
        "Target column `",
        target,
        "` was not found."
      )
    )
  }
  
  
  if (anyNA(data[[target]])) {
    
    stop(
      paste0(
        "Target column `",
        target,
        "` contains missing values."
      )
    )
  }
  
  
  method_control <- match.arg(
    method_control,
    c("cv", "oob")
  )
  
  
  if (
    method_control == "cv" &&
    (
      !is.numeric(cv_folds) ||
      length(cv_folds) != 1 ||
      cv_folds < 2 ||
      cv_folds %% 1 != 0
    )
  ) {
    
    stop(
      "`cv_folds` must be a single integer of at least 2."
    )
  }
  
  
  if (
    !is.numeric(n_repeats) ||
    length(n_repeats) != 1 ||
    n_repeats < 1 ||
    n_repeats %% 1 != 0
  ) {
    
    stop(
      "`n_repeats` must be a single positive integer."
    )
  }
  
  
  if (
    !is.numeric(n_cores) ||
    length(n_cores) != 1 ||
    n_cores < 1 ||
    n_cores %% 1 != 0
  ) {
    
    stop(
      "`n_cores` must be a single positive integer."
    )
  }
  
  
  if (
    !is.numeric(seed) ||
    length(seed) != 1 ||
    is.na(seed) ||
    seed < 0
  ) {
    
    stop(
      "`seed` must be a single non-negative integer."
    )
  }
  
  
  cv_folds <- as.integer(
    cv_folds
  )
  
  n_repeats <- as.integer(
    n_repeats
  )
  
  n_cores <- as.integer(
    n_cores
  )
  
  seed <- as.integer(
    seed
  )
  
  
  # ============================================================
  # 2. Copy input data
  # ============================================================
  
  # data.table can modify objects by reference.
  #
  # Work on a copy so AlgaeScan_train() can never alter the
  # user's original input object.
  
  data <- data.table::as.data.table(
    data.table::copy(data)
  )
  
  
  input_events <- nrow(
    data
  )
  
  
  input_class_counts <- table(
    as.character(
      data[[target]]
    )
  )
  
  
  # ============================================================
  # 3. Check downsampling arguments
  # ============================================================
  
  if (!is.null(downsample)) {
    
    if (
      !is.numeric(downsample) ||
      length(downsample) != 1 ||
      downsample <= 0 ||
      downsample %% 1 != 0
    ) {
      
      stop(
        "`downsample` must be a positive integer or NULL."
      )
    }
    
    
    downsample <- as.integer(
      downsample
    )
    
    
    if (
      is.null(downsample_by) ||
      length(downsample_by) == 0
    ) {
      
      stop(
        "`downsample_by` must be specified when downsampling is used."
      )
    }
    
    
    missing_groups <- downsample_by[
      !downsample_by %in% colnames(data)
    ]
    
    
    if (length(missing_groups) > 0) {
      
      stop(
        paste0(
          "Downsampling column(s) not found: ",
          paste(
            missing_groups,
            collapse = ", "
          )
        )
      )
    }
    
    
    for (column in downsample_by) {
      
      if (anyNA(data[[column]])) {
        
        stop(
          paste0(
            "Downsampling column `",
            column,
            "` contains missing values."
          )
        )
      }
    }
  }
  
  
  # ============================================================
  # 4. Determine predictor columns
  # ============================================================
  
  if (is.null(feature_cols)) {
    
    numeric_cols <- vapply(
      data,
      is.numeric,
      logical(1)
    )
    
    
    feature_cols <- colnames(
      data
    )[numeric_cols]
    
    
    feature_cols <- setdiff(
      feature_cols,
      target
    )
    
  } else {
    
    if (target %in% feature_cols) {
      
      stop(
        "`target` cannot also be included in `feature_cols`."
      )
    }
    
    
    missing_features <- feature_cols[
      !feature_cols %in% colnames(data)
    ]
    
    
    if (length(missing_features) > 0) {
      
      stop(
        paste0(
          "Predictor column(s) not found: ",
          paste(
            missing_features,
            collapse = ", "
          )
        )
      )
    }
    
    
    non_numeric <- feature_cols[
      !vapply(
        data[
          ,
          feature_cols,
          with = FALSE
        ],
        is.numeric,
        logical(1)
      )
    ]
    
    
    if (length(non_numeric) > 0) {
      
      stop(
        paste0(
          "Predictors must be numeric. Non-numeric column(s): ",
          paste(
            non_numeric,
            collapse = ", "
          )
        )
      )
    }
  }
  
  
  if (length(feature_cols) == 0) {
    
    stop(
      "No numeric predictor columns were found."
    )
  }
  
  
  # Random Forest-specific safety check.
  #
  # This is placed here because the number of selected predictors
  # is known only after feature selection.
  
  for (configuration in configurations) {
    
    if ("mtry" %in% colnames(configuration$tune_grid)) {
      
      invalid_mtry <- configuration$tune_grid$mtry[
        configuration$tune_grid$mtry >
          length(feature_cols)
      ]
      
      
      if (length(invalid_mtry) > 0) {
        
        stop(
          paste0(
            "`mtry` cannot be larger than the number of predictors (",
            length(feature_cols),
            "). Invalid value(s): ",
            paste(
              unique(invalid_mtry),
              collapse = ", "
            )
          )
        )
      }
    }
  }
  
  
  # ============================================================
  # 5. Prepare response labels
  # ============================================================
  
  original_classes <- unique(
    as.character(
      data[[target]]
    )
  )
  
  
  if (length(original_classes) < 2) {
    
    stop(
      "`target` must contain at least two classes."
    )
  }
  
  
  # caret requires syntactically valid class names for reliable
  # probability prediction.
  #
  # Example:
  #
  # non-algae -> non_algae
  
  safe_classes <- make.names(
    gsub(
      "-",
      "_",
      original_classes,
      fixed = TRUE
    ),
    unique = TRUE
  )
  
  
  original_to_safe <- stats::setNames(
    safe_classes,
    original_classes
  )
  
  
  safe_to_original <- stats::setNames(
    original_classes,
    safe_classes
  )
  
  
  # ============================================================
  # 6. Training summary
  # ============================================================
  
  message("")
  message("AlgaeScan model training")
  message("------------------------")
  
  message(
    "Method: ",
    method_name
  )
  
  message(
    "Target: ",
    target
  )
  
  message(
    "Input events: ",
    format(
      input_events,
      big.mark = ","
    )
  )
  
  message(
    "Predictor variables: ",
    length(feature_cols)
  )
  
  
  if (is.null(downsample)) {
    
    message(
      "Downsampling: none"
    )
    
  } else {
    
    message(
      "Downsampling: ",
      downsample,
      " per ",
      paste(
        downsample_by,
        collapse = " / "
      )
    )
  }
  
  
  if (method_control == "cv") {
    
    message(
      "Cross-validation: ",
      cv_folds,
      "-fold"
    )
    
  } else {
    
    message(
      "Validation: out-of-bag"
    )
  }
  
  
  message(
    "Training repetitions: ",
    n_repeats
  )
  
  message(
    "CPU cores: ",
    n_cores
  )
  
  
  if (is.null(preprocess)) {
    
    message(
      "Preprocessing: none"
    )
    
  } else {
    
    message(
      "Preprocessing: ",
      paste(
        preprocess,
        collapse = ", "
      )
    )
  }
  
  
  message("")
  
  
  # ============================================================
  # 7. Start parallel backend
  # ============================================================
  
  cluster <- NULL
  
  
  if (n_cores > 1) {
    
    cluster <- parallel::makePSOCKcluster(
      n_cores
    )
    
    
    doParallel::registerDoParallel(
      cluster
    )
    
  } else {
    
    foreach::registerDoSEQ()
  }
  
  
  # Always close the cluster, including when training fails.
  
  on.exit({
    
    if (!is.null(cluster)) {
      
      parallel::stopCluster(
        cluster
      )
    }
    
    
    foreach::registerDoSEQ()
    
  }, add = TRUE)
  
  
  # ============================================================
  # 8. Initialise training objects
  # ============================================================
  
  start_time <- Sys.time()
  
  
  best_model <- NULL
  
  best_accuracy <- -Inf
  
  best_parameters <- NULL
  
  best_repetition <- NULL
  
  best_training_events <- NULL
  
  best_training_class_counts <- NULL
  
  best_model_seed <- NULL
  
  
  all_results <- vector(
    mode = "list",
    length = n_repeats *
      length(configurations)
  )
  
  
  repetition_results <- vector(
    mode = "list",
    length = n_repeats
  )
  
  
  result_counter <- 0L
  
  
  # ============================================================
  # 9. Training repetitions
  # ============================================================
  
  for (repetition in seq_len(n_repeats)) {
    
    # ----------------------------------------------------------
    # Separate random seeds for independent stochastic processes
    # ----------------------------------------------------------
    
    repetition_seed <- as.integer(
      seed + repetition - 1L
    )
    
    
    downsample_seed <- repetition_seed
    
    
    cv_seed <- as.integer(
      repetition_seed + 1000L
    )
    
    
    model_seed <- as.integer(
      repetition_seed + 2000L
    )
    
    
    message(
      "========================================"
    )
    
    message(
      "Training repetition ",
      repetition,
      " of ",
      n_repeats
    )
    
    message(
      "========================================"
    )
    
    
    # ==========================================================
    # 9a. Optional downsampling
    # ==========================================================
    
    training_data <- data
    
    
    if (!is.null(downsample)) {
      
      message(
        "Downsampling training data..."
      )
      
      
      set.seed(
        downsample_seed
      )
      
      
      # --------------------------------------------------------
      # Generate grouping variable
      # --------------------------------------------------------
      
      if (length(downsample_by) == 1) {
        
        group_id <- training_data[[downsample_by]]
        
      } else {
        
        group_data <- training_data[
          ,
          downsample_by,
          with = FALSE
        ]
        
        
        group_id <- do.call(
          interaction,
          c(
            as.list(group_data),
            list(
              drop = TRUE
            )
          )
        )
        
        
        rm(
          group_data
        )
      }
      
      
      group_indices <- split(
        seq_len(
          nrow(training_data)
        ),
        group_id
      )
      
      
      # --------------------------------------------------------
      # Sample within each group
      # --------------------------------------------------------
      
      selected_rows <- unlist(
        lapply(
          group_indices,
          function(rows) {
            
            # If a group contains fewer events than requested,
            # retain all events without unnecessary randomization.
            
            if (length(rows) <= downsample) {
              
              return(
                rows
              )
            }
            
            
            sample(
              rows,
              size = downsample,
              replace = FALSE
            )
          }
        ),
        use.names = FALSE
      )
      
      
      # --------------------------------------------------------
      # Preserve original relative row order
      # --------------------------------------------------------
      
      selected_rows <- sort(
        selected_rows
      )
      
      
      training_data <- training_data[
        selected_rows,
      ]
      
      
      rm(
        group_id,
        group_indices,
        selected_rows
      )
    }
    
    
    n_training_events <- nrow(
      training_data
    )
    
    
    message(
      "Training events: ",
      format(
        n_training_events,
        big.mark = ","
      )
    )
    
    
    # ==========================================================
    # 9b. Extract predictors and response
    # ==========================================================
    
    X_train <- training_data[
      ,
      feature_cols,
      with = FALSE
    ]
    
    
    Y_original <- as.character(
      training_data[[target]]
    )
    
    
    training_class_counts <- table(
      Y_original
    )
    
    
    Y_train <- factor(
      original_to_safe[
        Y_original
      ],
      levels = safe_classes
    )
    
    
    # Every class must contain enough events to participate in
    # every cross-validation fold.
    
    if (
      method_control == "cv" &&
      any(
        table(Y_train) < cv_folds
      )
    ) {
      
      stop(
        paste0(
          "Each response class must contain at least ",
          cv_folds,
          " events for ",
          cv_folds,
          "-fold cross-validation."
        )
      )
    }
    
    
    rm(
      training_data,
      Y_original
    )
    
    
    # ----------------------------------------------------------
    # Report predictor memory size
    # ----------------------------------------------------------
    
    x_size_gb <- as.numeric(
      object.size(X_train)
    ) / 1024^3
    
    
    message(
      "Predictor data size: ",
      round(
        x_size_gb,
        2
      ),
      " GiB"
    )
    
    
    if (
      n_cores > 4 &&
      x_size_gb > 0.5
    ) {
      
      warning(
        paste0(
          "Large predictor dataset combined with ",
          n_cores,
          " parallel workers may require substantial memory."
        ),
        call. = FALSE
      )
    }
    
    
    # ==========================================================
    # 9c. Create explicit cross-validation folds
    # ==========================================================
    
    if (method_control == "cv") {
      
      set.seed(
        cv_seed
      )
      
      
      cv_index <- caret::createFolds(
        Y_train,
        k = cv_folds,
        returnTrain = TRUE
      )
      
    } else {
      
      cv_index <- NULL
    }
    
    
    # ==========================================================
    # 9d. Train each model configuration
    # ==========================================================
    
    repetition_best_accuracy <- -Inf
    
    repetition_best_parameters <- NULL
    
    
    for (
      config_index in seq_along(configurations)
    ) {
      
      configuration <- configurations[[config_index]]
      
      
      message("")
      
      message(
        "Configuration ",
        config_index,
        " of ",
        length(configurations),
        ": ",
        configuration$label
      )
      
      
      # --------------------------------------------------------
      # Build reproducible caret control
      # --------------------------------------------------------
      
      if (method_control == "cv") {
        
        caret_seeds <- .AlgaeScan_make_caret_seeds(
          seed = model_seed,
          n_resamples = cv_folds,
          n_models = nrow(
            configuration$tune_grid
          )
        )
        
        
        train_control <- caret::trainControl(
          method = "cv",
          number = cv_folds,
          index = cv_index,
          seeds = caret_seeds,
          trim = TRUE,
          returnData = FALSE,
          returnResamp = "none",
          savePredictions = FALSE,
          allowParallel = n_cores > 1
        )
        
      } else {
        
        caret_seeds <- NULL
        
        
        train_control <- caret::trainControl(
          method = "oob",
          trim = TRUE,
          returnData = FALSE,
          returnResamp = "none",
          savePredictions = FALSE,
          allowParallel = n_cores > 1
        )
      }
      
      
      # --------------------------------------------------------
      # Build caret::train() arguments
      # --------------------------------------------------------
      
      train_arguments <- list(
        x = X_train,
        y = Y_train,
        method = method,
        metric = "Accuracy",
        preProcess = preprocess,
        trControl = train_control,
        tuneGrid = configuration$tune_grid
      )
      
      
      # Add model-specific arguments such as:
      #
      # RF:
      #   ntree
      #
      # NNET:
      #   trace
      #   maxit
      #   MaxNWts
      
      train_arguments <- c(
        train_arguments,
        configuration$train_args
      )
      
      
      # --------------------------------------------------------
      # Train configuration
      # --------------------------------------------------------
      
      set.seed(
        model_seed
      )
      
      
      model_i <- do.call(
        caret::train,
        train_arguments
      )
      
      
      # ========================================================
      # 9e. Recover selected tuning parameters
      # ========================================================
      
      best_tune <- as.list(
        model_i$bestTune[
          1,
          ,
          drop = FALSE
        ]
      )
      
      
      # --------------------------------------------------------
      # Rename caret-specific parameter names when needed
      #
      # Example:
      # size -> hidden_units
      # --------------------------------------------------------
      
      if (!is.null(parameter_map)) {
        
        for (
          old_name in names(parameter_map)
        ) {
          
          if (old_name %in% names(best_tune)) {
            
            names(best_tune)[
              names(best_tune) == old_name
            ] <- parameter_map[[old_name]]
          }
        }
      }
      
      
      # Combine fixed parameters and tuned parameters.
      #
      # Example RF:
      #
      # n_tree = 10
      # mtry = 50
      
      parameters_i <- c(
        configuration$fixed_parameters,
        best_tune
      )
      
      
      # ========================================================
      # 9f. Store validation results
      # ========================================================
      
      results_i <- data.table::as.data.table(
        data.table::copy(
          model_i$results
        )
      )
      
      
      # Rename caret-specific tuning parameter columns.
      
      if (!is.null(parameter_map)) {
        
        for (
          old_name in names(parameter_map)
        ) {
          
          if (old_name %in% colnames(results_i)) {
            
            data.table::setnames(
              results_i,
              old_name,
              parameter_map[[old_name]]
            )
          }
        }
      }
      
      
      # Add fixed model parameters such as n_tree.
      
      if (
        length(
          configuration$fixed_parameters
        ) > 0
      ) {
        
        for (
          parameter_name in
          names(
            configuration$fixed_parameters
          )
        ) {
          
          results_i[
            ,
            (parameter_name) :=
              configuration$fixed_parameters[[parameter_name]]
          ]
        }
      }
      
      
      # Add reproducibility information.
      
      results_i[
        ,
        `:=`(
          method = method,
          repetition = repetition,
          repetition_seed = repetition_seed,
          downsample_seed = downsample_seed,
          cv_seed = if (
            method_control == "cv"
          ) cv_seed else NA_integer_,
          model_seed = model_seed
        )
      ]
      
      
      result_counter <- result_counter + 1L
      
      
      all_results[[result_counter]] <- results_i
      
      
      # ========================================================
      # 9g. Compare model performance
      # ========================================================
      
      accuracy_i <- max(
        model_i$results$Accuracy,
        na.rm = TRUE
      )
      
      
      if (!is.finite(accuracy_i)) {
        
        stop(
          paste0(
            "Training failed to generate a valid accuracy for ",
            configuration$label,
            "."
          )
        )
      }
      
      
      message(
        "Best parameters: ",
        paste(
          paste0(
            names(parameters_i),
            " = ",
            unlist(parameters_i)
          ),
          collapse = ", "
        )
      )
      
      
      message(
        "Best accuracy: ",
        round(
          accuracy_i,
          6
        )
      )
      
      
      # --------------------------------------------------------
      # Best model within this repetition
      # --------------------------------------------------------
      
      if (
        accuracy_i >
        repetition_best_accuracy
      ) {
        
        repetition_best_accuracy <-
          accuracy_i
        
        repetition_best_parameters <-
          parameters_i
      }
      
      
      # --------------------------------------------------------
      # Best model across all repetitions
      #
      # If two models have exactly the same accuracy, the first
      # encountered model is retained.
      # --------------------------------------------------------
      
      if (
        accuracy_i >
        best_accuracy
      ) {
        
        best_accuracy <- accuracy_i
        
        best_parameters <- parameters_i
        
        best_repetition <- repetition
        
        best_training_events <-
          n_training_events
        
        best_training_class_counts <-
          training_class_counts
        
        best_model_seed <- model_seed
        
        
        # Replace previously stored best model.
        
        if (!is.null(best_model)) {
          
          rm(
            best_model
          )
        }
        
        
        best_model <- model_i
        
        
        rm(
          model_i
        )
        
      } else {
        
        # Do not retain full caret objects for models that were
        # not selected.
        
        rm(
          model_i
        )
      }
    }
    
    
    # ==========================================================
    # 9h. Store repetition summary
    # ==========================================================
    
    repetition_row <- data.table::data.table(
      repetition = repetition,
      repetition_seed = repetition_seed,
      downsample_seed = downsample_seed,
      cv_seed = if (
        method_control == "cv"
      ) cv_seed else NA_integer_,
      model_seed = model_seed,
      training_events = n_training_events,
      best_accuracy = repetition_best_accuracy
    )
    
    
    # Add best hyperparameters from this repetition.
    
    for (
      parameter_name in
      names(
        repetition_best_parameters
      )
    ) {
      
      output_name <- paste0(
        "best_",
        parameter_name
      )
      
      
      repetition_row[
        ,
        (output_name) :=
          repetition_best_parameters[[parameter_name]]
      ]
    }
    
    
    repetition_results[[repetition]] <- repetition_row
    
    
    message("")
    
    message(
      "Repetition ",
      repetition,
      " complete."
    )
    
    
    message(
      "Best accuracy: ",
      round(
        repetition_best_accuracy,
        6
      )
    )
    
    
    # ----------------------------------------------------------
    # Remove large repetition-specific objects
    # ----------------------------------------------------------
    
    rm(
      X_train,
      Y_train,
      training_class_counts,
      cv_index
    )
    
    
    gc(
      verbose = FALSE
    )
    
    
    message("")
  }
  
  
  # ============================================================
  # 10. Combine validation results
  # ============================================================
  
  validation_results <- data.table::rbindlist(
    all_results[
      seq_len(
        result_counter
      )
    ],
    use.names = TRUE,
    fill = TRUE
  )
  
  
  repetition_results <- data.table::rbindlist(
    repetition_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  
  end_time <- Sys.time()
  
  
  # ============================================================
  # 11. Record software versions
  # ============================================================
  
  software_versions <- list(
    
    R = R.version.string,
    
    AlgaeScan =
      .AlgaeScan_package_version(
        "AlgaeScan"
      ),
    
    caret =
      .AlgaeScan_package_version(
        "caret"
      ),
    
    data.table =
      .AlgaeScan_package_version(
        "data.table"
      ),
    
    foreach =
      .AlgaeScan_package_version(
        "foreach"
      ),
    
    doParallel =
      .AlgaeScan_package_version(
        "doParallel"
      )
  )
  
  
  # Add the model-specific backend package.
  
  if (!is.null(backend_package)) {
    
    software_versions[[backend_package]] <- .AlgaeScan_package_version(backend_package)
  }
  
  
  # ============================================================
  # 12. Build AlgaeScan model object
  # ============================================================
  
  output <- list(
    
    # Fitted caret model
    model = best_model,
    
    # Model information
    method = method,
    target = target,
    features = feature_cols,
    
    # Original biological class names
    classes = original_classes,
    
    # Mapping from internal caret labels back to original labels
    label_map = safe_to_original,
    
    # Best hyperparameters
    best_parameters = best_parameters,
    
    # Validation performance
    best_accuracy = best_accuracy,
    
    # Repetition selected
    best_repetition = best_repetition,
    
    # Full validation results
    validation_results = validation_results,
    
    # Summary of each repetition
    repetition_results = repetition_results,
    
    # Training dataset information
    input_events = input_events,
    input_class_counts = input_class_counts,
    
    training_events = best_training_events,
    training_class_counts =
      best_training_class_counts,
    
    # Training settings
    settings = list(
      method = method,
      tuning = tuning_settings,
      downsample = downsample,
      downsample_by = downsample_by,
      cv_folds = cv_folds,
      n_repeats = n_repeats,
      method_control = method_control,
      preprocess = preprocess,
      n_cores = n_cores,
      seed = seed,
      preserve_input_order = TRUE,
      explicit_cv_folds =
        method_control == "cv",
      explicit_caret_seeds =
        method_control == "cv"
    ),
    
    # Reproducibility metadata
    reproducibility = list(
      
      best_repetition_seed =
        as.integer(
          seed +
            best_repetition -
            1L
        ),
      
      best_downsample_seed =
        as.integer(
          seed +
            best_repetition -
            1L
        ),
      
      best_cv_seed =
        if (
          method_control == "cv"
        ) {
          as.integer(
            seed +
              best_repetition -
              1L +
              1000L
          )
        } else {
          NA_integer_
        },
      
      best_model_seed =
        best_model_seed,
      
      software_versions =
        software_versions
    ),
    
    # Total elapsed training time
    training_time =
      end_time -
      start_time
  )
  
  
  class(
    output
  ) <- "AlgaeScan_model"
  
  
  # ============================================================
  # 13. Training summary
  # ============================================================
  
  message(
    "AlgaeScan training complete."
  )
  
  
  message(
    "Method: ",
    method_name
  )
  
  
  message(
    "Selected repetition: ",
    best_repetition,
    " of ",
    n_repeats
  )
  
  
  message(
    "Selected parameters: ",
    paste(
      paste0(
        names(best_parameters),
        " = ",
        unlist(best_parameters)
      ),
      collapse = ", "
    )
  )
  
  
  message(
    "Best validation accuracy: ",
    round(
      best_accuracy,
      6
    )
  )
  
  
  return(
    output
  )
}

#' Shared unsupervised training engine for AlgaeScan
#'
#' Internal training engine shared by the unsupervised and one-class
#' AlgaeScan training functions.
#'
#' The function handles functionality shared across unsupervised
#' algorithms, including:
#'
#' \itemize{
#'   \item predictor selection
#'   \item optional downsampling
#'   \item repeated training
#'   \item cross-validation folds
#'   \item optional preprocessing
#'   \item model-specific fitting and prediction
#'   \item validation scoring
#'   \item model selection
#'   \item reproducibility metadata
#' }
#'
#' One-Class SVM and Isolation Forest are trained using only events
#' belonging to the reference class. Validation is performed using
#' all classes.
#'
#' Gaussian Mixture Models are trained using all events without using
#' the class labels. Reference labels are used only for validation.
#'
#' Isolation Forest fitting and prediction are performed using the
#' Python implementation distributed with AlgaeScan.
#'
#' @param data Training data as a `data.frame` or `data.table`.
#'
#' @param target Character string giving the column containing the
#'   reference labels used for validation.
#'
#' @param method Character string identifying the unsupervised model.
#'   Supported values are `"one_class_svm"`, `"gmm"` and
#'   `"isoforest"`.
#'
#' @param method_name Human-readable model name.
#'
#' @param feature_cols Optional character vector containing predictor
#'   columns. If `NULL`, numeric predictor columns are selected
#'   automatically.
#'
#' @param downsample Optional integer giving the maximum number of
#'   events retained per downsampling group. If `NULL`, all events
#'   are used.
#'
#' @param downsample_by Character vector containing the column or
#'   columns defining downsampling groups.
#'
#' @param cv_folds Number of cross-validation folds.
#'
#' @param n_repeats Number of complete training repetitions.
#'
#' @param preprocess Optional preprocessing applied before
#'   cross-validation. Typically `c("center", "scale")`.
#'
#' @param n_cores Number of CPU cores requested through the common
#'   AlgaeScan training interface.
#'
#' @param seed Integer controlling reproducible downsampling,
#'   cross-validation and model fitting.
#'
#' @param ref_label Optional character string defining the reference
#'   or normal class for one-class methods.
#'
#' @param fold_parameter_name Optional character string giving the
#'   name of a model parameter whose values are assigned sequentially
#'   across cross-validation folds.
#'
#' @param fold_parameter_values Optional vector containing one model
#'   parameter value or one value per cross-validation fold.
#'
#' @param score_metric Character string describing the validation
#'   metric.
#'
#' @param python_random_state Integer passed to the Python Isolation
#'   Forest implementation as `random_state`.
#'
#' @return An object of class `AlgaeScan_model` containing the selected
#'   fitted model, validation results, repetition summaries, selected
#'   parameters, preprocessing information and reproducibility
#'   metadata.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_train_unsupervised <- function(
    data,
    target,
    method,
    method_name,
    feature_cols = NULL,
    downsample = NULL,
    downsample_by = "Species",
    cv_folds = 5,
    n_repeats = 1,
    preprocess = c("center", "scale"),
    n_cores = 1,
    seed = 123,
    ref_label = NULL,
    fold_parameter_name = NULL,
    fold_parameter_values = NULL,
    score_metric = "Accuracy",
    python_random_state = 42
) {
  
  # ============================================================
  # 1. General input checks
  # ============================================================
  
  if (!is.data.frame(data)) {
    
    stop(
      "`data` must be a data.frame or data.table."
    )
  }
  
  
  if (nrow(data) == 0) {
    
    stop(
      "`data` contains no events."
    )
  }
  
  
  if (
    !is.character(target) ||
    length(target) != 1 ||
    is.na(target) ||
    target == ""
  ) {
    
    stop(
      "`target` must be a single column name."
    )
  }
  
  
  if (!target %in% colnames(data)) {
    
    stop(
      paste0(
        "Target column `",
        target,
        "` was not found."
      )
    )
  }
  
  
  if (anyNA(data[[target]])) {
    
    stop(
      paste0(
        "Target column `",
        target,
        "` contains missing values."
      )
    )
  }
  
  
  method <- match.arg(
    method,
    c(
      "one_class_svm",
      "gmm",
      "isoforest"
    )
  )
  
  
  if (
    !is.numeric(cv_folds) ||
    length(cv_folds) != 1 ||
    cv_folds < 2 ||
    cv_folds %% 1 != 0
  ) {
    
    stop(
      "`cv_folds` must be a single integer of at least 2."
    )
  }
  
  
  if (
    !is.numeric(n_repeats) ||
    length(n_repeats) != 1 ||
    n_repeats < 1 ||
    n_repeats %% 1 != 0
  ) {
    
    stop(
      "`n_repeats` must be a single positive integer."
    )
  }
  
  
  if (
    !is.numeric(n_cores) ||
    length(n_cores) != 1 ||
    n_cores < 1 ||
    n_cores %% 1 != 0
  ) {
    
    stop(
      "`n_cores` must be a single positive integer."
    )
  }
  
  
  if (!is.null(downsample)) {
    
    if (
      !is.numeric(downsample) ||
      length(downsample) != 1 ||
      downsample < 1 ||
      downsample %% 1 != 0
    ) {
      
      stop(
        "`downsample` must be a positive integer or NULL."
      )
    }
    
    
    missing_downsample_cols <- downsample_by[
      !downsample_by %in% colnames(data)
    ]
    
    
    if (length(missing_downsample_cols) > 0) {
      
      stop(
        paste0(
          "Downsampling column(s) not found: ",
          paste(
            missing_downsample_cols,
            collapse = ", "
          )
        )
      )
    }
  }
  
  
  # ============================================================
  # 2. Copy data
  # ============================================================
  
  data <- data.table::as.data.table(
    data.table::copy(data)
  )
  
  
  input_events <- nrow(
    data
  )
  
  
  input_class_counts <- table(
    data[[target]]
  )
  
  
  # ============================================================
  # 3. Determine predictor columns
  # ============================================================
  
  if (is.null(feature_cols)) {
    
    numeric_cols <- vapply(
      data,
      is.numeric,
      logical(1)
    )
    
    
    feature_cols <- colnames(
      data
    )[numeric_cols]
    
    
    feature_cols <- setdiff(
      feature_cols,
      target
    )
    
  } else {
    
    if (target %in% feature_cols) {
      
      stop(
        "`target` cannot also be included in `feature_cols`."
      )
    }
    
    
    missing_features <- feature_cols[
      !feature_cols %in% colnames(data)
    ]
    
    
    if (length(missing_features) > 0) {
      
      stop(
        paste0(
          "Predictor column(s) not found: ",
          paste(
            missing_features,
            collapse = ", "
          )
        )
      )
    }
    
    
    non_numeric <- feature_cols[
      !vapply(
        data[
          ,
          feature_cols,
          with = FALSE
        ],
        is.numeric,
        logical(1)
      )
    ]
    
    
    if (length(non_numeric) > 0) {
      
      stop(
        paste0(
          "Predictors must be numeric. Non-numeric column(s): ",
          paste(
            non_numeric,
            collapse = ", "
          )
        )
      )
    }
  }
  
  
  if (length(feature_cols) == 0) {
    
    stop(
      "No numeric predictor columns were found."
    )
  }
  
  
  # ============================================================
  # 4. Check reference classes
  # ============================================================
  
  classes <- unique(
    as.character(
      data[[target]]
    )
  )
  
  
  if (length(classes) < 2) {
    
    stop(
      "`target` must contain at least two classes."
    )
  }
  
  
  other_label <- NULL
  
  
  if (
    method %in%
    c(
      "one_class_svm",
      "isoforest"
    )
  ) {
    
    if (
      is.null(ref_label) ||
      !ref_label %in% classes
    ) {
      
      stop(
        paste0(
          "`ref_label` must identify a class present in `",
          target,
          "`."
        )
      )
    }
    
    
    other_classes <- setdiff(
      classes,
      ref_label
    )
    
    
    if (length(other_classes) != 1) {
      
      stop(
        paste0(
          method_name,
          " currently requires exactly two reference classes."
        )
      )
    }
    
    
    other_label <- other_classes[1]
  }
  
  
  # ============================================================
  # 5. Check fold-specific parameters
  # ============================================================
  
  if (!is.null(fold_parameter_values)) {
    
    if (length(fold_parameter_values) == 1) {
      
      fold_parameter_values <- rep(
        fold_parameter_values,
        cv_folds
      )
    }
    
    
    if (
      length(fold_parameter_values) !=
      cv_folds
    ) {
      
      stop(
        paste0(
          "`",
          fold_parameter_name,
          "` must contain either one value or one value per ",
          "cross-validation fold."
        )
      )
    }
  }
  
  
  # ============================================================
  # 6. Load Python Isolation Forest implementation
  # ============================================================
  
  isoforest_python <- NULL
  
  
  if (method == "isoforest") {
    
    isoforest_python <-
      .AlgaeScan_load_isoforest_python()
  }
  
  
  # ============================================================
  # 7. Training summary
  # ============================================================
  
  message("")
  message("AlgaeScan model training")
  message("------------------------")
  
  message(
    "Method: ",
    method_name
  )
  
  message(
    "Target: ",
    target
  )
  
  message(
    "Input events: ",
    format(
      input_events,
      big.mark = ","
    )
  )
  
  message(
    "Predictor variables: ",
    length(feature_cols)
  )
  
  
  if (is.null(downsample)) {
    
    message(
      "Downsampling: none"
    )
    
  } else {
    
    message(
      "Downsampling: ",
      downsample,
      " per ",
      paste(
        downsample_by,
        collapse = " / "
      )
    )
  }
  
  
  message(
    "Cross-validation: ",
    cv_folds,
    "-fold"
  )
  
  message(
    "Training repetitions: ",
    n_repeats
  )
  
  
  if (is.null(preprocess)) {
    
    message(
      "Preprocessing: none"
    )
    
  } else {
    
    message(
      "Preprocessing: ",
      paste(
        preprocess,
        collapse = ", "
      )
    )
  }
  
  
  if (!is.null(ref_label)) {
    
    message(
      "Reference class: ",
      ref_label
    )
  }
  
  
  message("")
  
  
  # ============================================================
  # 8. Initialise result containers
  # ============================================================
  
  start_time <- Sys.time()
  
  
  all_results <- list()
  
  repetition_results <- list()
  
  
  result_counter <- 0L
  
  
  best_model <- NULL
  
  best_score <- -Inf
  
  best_parameters <- list()
  
  best_repetition <- NULL
  
  best_fold <- NULL
  
  best_preprocess <- NULL
  
  best_training_events <- NULL
  
  best_training_class_counts <- NULL
  
  
  # ============================================================
  # 9. Training repetitions
  # ============================================================
  
  for (repetition in seq_len(n_repeats)) {
    
    repetition_seed <- as.integer(
      seed + repetition - 1L
    )
    
    
    downsample_seed <- repetition_seed
    
    
    cv_seed <- as.integer(
      repetition_seed + 1000L
    )
    
    
    model_seed <- as.integer(
      repetition_seed + 2000L
    )
    
    
    message(
      "========================================"
    )
    
    message(
      "Repetition ",
      repetition,
      " of ",
      n_repeats
    )
    
    message(
      "========================================"
    )
    
    
    # ==========================================================
    # 9a. Optional downsampling
    # ==========================================================
    
    training_data <- data
    
    
    if (!is.null(downsample)) {
      
      message(
        "Downsampling training data..."
      )
      
      
      set.seed(
        downsample_seed
      )
      
      
      if (length(downsample_by) == 1) {
        
        group_id <-
          training_data[[downsample_by]]
        
      } else {
        
        group_data <- training_data[
          ,
          downsample_by,
          with = FALSE
        ]
        
        
        group_id <- do.call(
          interaction,
          c(
            as.list(group_data),
            list(
              drop = TRUE
            )
          )
        )
      }
      
      
      group_indices <- split(
        seq_len(
          nrow(training_data)
        ),
        group_id
      )
      
      
      selected_rows <- unlist(
        lapply(
          group_indices,
          function(rows) {
            
            if (length(rows) <= downsample) {
              
              return(
                rows
              )
            }
            
            
            sample(
              rows,
              size = downsample,
              replace = FALSE
            )
          }
        ),
        use.names = FALSE
      )
      
      
      # Preserve original relative row order.
      
      selected_rows <- sort(
        selected_rows
      )
      
      
      training_data <-
        training_data[
          selected_rows,
        ]
    }
    
    
    n_training_events <- nrow(
      training_data
    )
    
    
    training_class_counts <- table(
      training_data[[target]]
    )
    
    
    message(
      "Training events: ",
      format(
        n_training_events,
        big.mark = ","
      )
    )
    
    
    # ==========================================================
    # 9b. Extract predictors and response
    # ==========================================================
    
    X_train <- as.data.frame(
      training_data[
        ,
        feature_cols,
        with = FALSE
      ]
    )
    
    
    Y_train <- as.character(
      training_data[[target]]
    )
    
    
    rm(
      training_data
    )
    
    
    # ==========================================================
    # 9c. Preprocess predictors
    # ==========================================================
    
    preprocess_model <- NULL
    
    
    if (!is.null(preprocess)) {
      
      preprocess_model <- caret::preProcess(
        X_train,
        method = preprocess
      )
      
      
      X_train <- stats::predict(
        preprocess_model,
        X_train
      )
    }
    
    
    # ==========================================================
    # 9d. Create cross-validation folds
    # ==========================================================
    
    set.seed(
      cv_seed
    )
    
    
    cv_index <- caret::createFolds(
      Y_train,
      k = cv_folds,
      returnTrain = TRUE
    )
    
    
    repetition_best_score <- -Inf
    
    repetition_best_fold <- NULL
    
    repetition_best_parameter <- NULL
    
    
    # ==========================================================
    # 9e. Train each fold
    # ==========================================================
    
    for (fold in seq_len(cv_folds)) {
      
      message(
        "Fold ",
        fold,
        " of ",
        cv_folds
      )
      
      
      train_index <-
        cv_index[[fold]]
      
      
      validation_index <- setdiff(
        seq_len(
          nrow(X_train)
        ),
        train_index
      )
      
      
      X_fold_train <- X_train[
        train_index,
        ,
        drop = FALSE
      ]
      
      
      X_fold_validation <- X_train[
        validation_index,
        ,
        drop = FALSE
      ]
      
      
      Y_fold_train <-
        Y_train[train_index]
      
      
      Y_fold_validation <-
        Y_train[validation_index]
      
      
      parameter_value <- NULL
      
      
      if (!is.null(fold_parameter_values)) {
        
        parameter_value <-
          fold_parameter_values[fold]
      }
      
      
      set.seed(
        model_seed + fold - 1L
      )
      
      
      # --------------------------------------------------------
      # One-Class SVM
      # --------------------------------------------------------
      
      if (method == "one_class_svm") {
        
        reference_rows <-
          Y_fold_train == ref_label
        
        
        model_fold <- e1071::svm(
          X_fold_train[
            reference_rows,
            ,
            drop = FALSE
          ],
          type = "one-classification",
          kernel = "radial",
          nu = parameter_value
        )
        
        
        prediction_raw <- stats::predict(
          model_fold,
          newdata = X_fold_validation
        )
        
        
        prediction <- ifelse(
          as.logical(prediction_raw),
          ref_label,
          other_label
        )
        
        
        score <- mean(
          prediction ==
            Y_fold_validation
        )
      }
      
      
      # --------------------------------------------------------
      # Gaussian Mixture Model
      # --------------------------------------------------------
      
      if (method == "gmm") {
        
        model_fold <- mclust::Mclust(
          X_fold_train
        )
        
        
        prediction <- stats::predict(
          model_fold,
          newdata = X_fold_validation
        )$classification
        
        
        score <- mclust::adjustedRandIndex(
          prediction,
          Y_fold_validation
        )
      }
      
      
      # --------------------------------------------------------
      # Isolation Forest
      # --------------------------------------------------------
      
      if (method == "isoforest") {
        
        reference_rows <-
          Y_fold_train == ref_label
        
        
        X_reference <- X_fold_train[
          reference_rows,
          ,
          drop = FALSE
        ]
        
        
        # Actual Isolation Forest training occurs in Python.
        
        model_fold <-
          isoforest_python$iso_forest_train(
            Xtrain = as.matrix(
              X_reference
            ),
            cont_val = parameter_value,
            random_state = as.integer(
              python_random_state
            )
          )
        
        
        # Prediction also occurs in Python.
        
        prediction_raw <-
          isoforest_python$iso_forest_predict(
            isof_model = model_fold,
            test_df = as.matrix(
              X_fold_validation
            )
          )
        
        
        # scikit-learn IsolationForest returns:
        #
        #  1 = inlier
        # -1 = outlier
        
        prediction <- ifelse(
          as.integer(prediction_raw) == 1L,
          ref_label,
          other_label
        )
        
        
        score <- mean(
          prediction ==
            Y_fold_validation
        )
      }
      
      
      # ========================================================
      # 9f. Store fold results
      # ========================================================
      
      result_row <- data.table::data.table(
        method = method,
        repetition = repetition,
        fold = fold,
        score = score,
        score_metric = score_metric,
        repetition_seed = repetition_seed,
        downsample_seed = downsample_seed,
        cv_seed = cv_seed,
        model_seed = model_seed
      )
      
      
      if (!is.null(fold_parameter_name)) {
        
        result_row[
          ,
          (fold_parameter_name) :=
            parameter_value
        ]
      }
      
      
      result_counter <-
        result_counter + 1L
      
      
      all_results[[result_counter]] <-
        result_row
      
      
      # ========================================================
      # 9g. Best model within repetition
      # ========================================================
      
      if (score > repetition_best_score) {
        
        repetition_best_score <- score
        
        repetition_best_fold <- fold
        
        repetition_best_parameter <-
          parameter_value
      }
      
      
      # ========================================================
      # 9h. Best model overall
      # ========================================================
      
      if (score > best_score) {
        
        best_score <- score
        
        best_model <- model_fold
        
        best_repetition <- repetition
        
        best_fold <- fold
        
        best_preprocess <-
          preprocess_model
        
        best_training_events <-
          n_training_events
        
        best_training_class_counts <-
          training_class_counts
        
        
        if (!is.null(fold_parameter_name)) {
          
          best_parameters <- stats::setNames(
            list(
              parameter_value
            ),
            fold_parameter_name
          )
          
        } else {
          
          best_parameters <- list()
        }
      }
    }
    
    
    # ==========================================================
    # 9i. Store repetition result
    # ==========================================================
    
    repetition_row <- data.table::data.table(
      repetition = repetition,
      repetition_seed = repetition_seed,
      downsample_seed = downsample_seed,
      cv_seed = cv_seed,
      model_seed = model_seed,
      training_events = n_training_events,
      best_score = repetition_best_score,
      best_fold = repetition_best_fold,
      score_metric = score_metric
    )
    
    
    if (!is.null(fold_parameter_name)) {
      
      repetition_row[
        ,
        (fold_parameter_name) :=
          repetition_best_parameter
      ]
    }
    
    
    repetition_results[[repetition]] <-
      repetition_row
  }
  
  
  # ============================================================
  # 10. Combine results
  # ============================================================
  
  validation_results <- data.table::rbindlist(
    all_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  
  repetition_results <- data.table::rbindlist(
    repetition_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  
  # ============================================================
  # 11. Build model object
  # ============================================================
  
  end_time <- Sys.time()
  
  
  output <- list(
    
    # Selected fitted model
    model = best_model,
    
    # Model information
    method = method,
    target = target,
    features = feature_cols,
    classes = classes,
    ref_label = ref_label,
    
    # Preprocessing
    preprocess_model = best_preprocess,
    
    # Selected parameters and validation score
    best_parameters = best_parameters,
    best_score = best_score,
    score_metric = score_metric,
    best_repetition = best_repetition,
    best_fold = best_fold,
    
    # Validation results
    validation_results = validation_results,
    repetition_results = repetition_results,
    
    # Input data summary
    input_events = input_events,
    input_class_counts = input_class_counts,
    
    # Selected training subset summary
    training_events = best_training_events,
    training_class_counts = best_training_class_counts,
    
    # Training settings
    settings = list(
      downsample = downsample,
      downsample_by = downsample_by,
      cv_folds = cv_folds,
      n_repeats = n_repeats,
      preprocess = preprocess,
      n_cores = n_cores,
      seed = seed,
      python_random_state = if (
        method == "isoforest"
      ) {
        python_random_state
      } else {
        NULL
      }
    ),
    
    # Timing
    training_time = as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    )
  )
  
  
  class(output) <- c(
    "AlgaeScan_model",
    "list"
  )
  
  
  return(
    output
  )
}

#' Generate deterministic caret seeds
#'
#' Internal helper used to create the seed structure expected by
#' `caret::trainControl()`.
#'
#' One seed vector is generated for each cross-validation resample,
#' followed by one additional seed used when fitting the final model
#' to all training events.
#'
#' @param seed Integer seed.
#' @param n_resamples Number of cross-validation resamples.
#' @param n_models Number of tuning combinations evaluated within
#'   each resample.
#'
#' @return A list of integer seed vectors suitable for the `seeds`
#'   argument of `caret::trainControl()`.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_make_caret_seeds <- function(
    seed,
    n_resamples,
    n_models
) {
  
  set.seed(
    seed
  )
  
  
  seeds <- vector(
    mode = "list",
    length = n_resamples + 1L
  )
  
  
  # One vector of model seeds per resample.
  
  for (
    i in seq_len(n_resamples)
  ) {
    
    seeds[[i]] <- sample.int(.Machine$integer.max,n_models)
  }
  
  
  # One additional seed for the final model fitted to the
  # complete training dataset.
  
  seeds[[n_resamples + 1L]] <- sample.int(
    .Machine$integer.max,
    1L
  )
  
  
  return(
    seeds
  )
}



#' Safely retrieve an R package version
#'
#' Internal helper used to record software versions in trained
#' AlgaeScan model objects.
#'
#' When a package version cannot be retrieved, `NA_character_` is
#' returned instead of interrupting model training.
#'
#' @param package Character string containing a package name.
#'
#' @return Character string containing the installed package version,
#'   or `NA_character_` if the version cannot be retrieved.
#'
#' @keywords internal
#' @noRd
.AlgaeScan_package_version <- function(
    package
) {
  
  tryCatch(
    
    as.character(
      utils::packageVersion(
        package
      )
    ),
    
    error = function(e) {
      
      NA_character_
    }
  )
}



#' Predict using an AlgaeScan model
#'
#' Applies a trained AlgaeScan model to new spectral flow cytometry
#' data.
#'
#' Predictor variables are automatically selected using the features
#' stored in the model and reordered to exactly match the order used
#' during training. Additional metadata columns in `new_data` are
#' ignored.
#'
#' Prediction behaviour depends on the fitted model:
#'
#' \itemize{
#'   \item Supervised `caret` models return predicted classes and,
#'   when requested, class probabilities.
#'
#'   \item One-Class SVM models return the stored reference class for
#'   inliers and the alternative class for outliers. When available,
#'   the SVM decision value is returned as `normality_score`.
#'
#'   \item Gaussian Mixture Models return mixture-component
#'   assignments rather than biological class labels. When requested,
#'   posterior mixture-component probabilities are also returned.
#'
#'   \item Isolation Forest models return the stored reference class
#'   for inliers and the alternative class for outliers. When
#'   requested, a normalized anomaly-detection score is returned as
#'   `normality_score`, where larger values indicate events that are
#'   more similar to the reference class.
#' }
#'
#' For supervised `caret` models, preprocessing learned during
#' training is applied automatically by the fitted `caret` model.
#'
#' For One-Class SVM and Isolation Forest models, preprocessing stored
#' during training is applied automatically before prediction.
#'
#' @param model An object of class `AlgaeScan_model`, generated by
#'   [AlgaeScan_train()] or [AlgaeScan_load_model()].
#'
#' @param new_data A `data.frame` or `data.table` containing the data
#'   to classify.
#'
#' @param probabilities Logical. If `TRUE`, additional model-specific
#'   probability or score information is returned when available.
#'   Default is `TRUE`.
#'
#' @param keep_data Logical. If `TRUE`, predictions are appended to a
#'   copy of `new_data`. If `FALSE`, only prediction results are
#'   returned. Default is `FALSE`.
#'
#' @param prediction_name Optional character string giving the name of
#'   the prediction column. If `NULL`, the name is generated from the
#'   model target. For GMM models the default is
#'   `"cluster_predicted"`.
#'
#' @param seed Integer used to make prediction reproducible.
#'   Default is `123`.
#'
#' @return A `data.table` containing predictions and, when requested
#'   and supported by the model, probabilities or model-specific
#'   scores. If `keep_data = TRUE`, these columns are appended to a
#'   copy of the input data.
#'
#' @examples
#' \dontrun{
#'
#' pred <- AlgaeScan_predict(
#'   model = model_A,
#'   new_data = prepared$test
#' )
#'
#' pred <- AlgaeScan_predict(
#'   model = model_A,
#'   new_data = prepared$test,
#'   keep_data = TRUE
#' )
#'
#' }
#'
#' @export
AlgaeScan_predict <- function(
    model,
    new_data,
    probabilities = TRUE,
    keep_data = FALSE,
    prediction_name = NULL,
    seed = 123
) {
  
  # ============================================================
  # 1. Check inputs
  # ============================================================
  
  if (!inherits(model, "AlgaeScan_model")) {
    
    stop(
      paste0(
        "`model` must be an `AlgaeScan_model` object.\n",
        "Use `AlgaeScan_load_model()` to import legacy caret models."
      )
    )
  }
  
  
  if (!is.data.frame(new_data)) {
    
    stop(
      "`new_data` must be a data.frame or data.table."
    )
  }
  
  
  if (nrow(new_data) == 0) {
    
    stop(
      "`new_data` contains no events."
    )
  }
  
  
  if (is.null(model$model)) {
    
    stop(
      "The AlgaeScan model does not contain a fitted model."
    )
  }
  
  
  if (
    !is.logical(probabilities) ||
    length(probabilities) != 1 ||
    is.na(probabilities)
  ) {
    
    stop(
      "`probabilities` must be TRUE or FALSE."
    )
  }
  
  
  if (
    !is.logical(keep_data) ||
    length(keep_data) != 1 ||
    is.na(keep_data)
  ) {
    
    stop(
      "`keep_data` must be TRUE or FALSE."
    )
  }
  
  
  if (
    !is.numeric(seed) ||
    length(seed) != 1 ||
    is.na(seed)
  ) {
    
    stop(
      "`seed` must be a single numeric value."
    )
  }
  
  
  new_data <- data.table::as.data.table(
    new_data
  )
  
  
  # ============================================================
  # 2. Select predictor variables
  # ============================================================
  
  features <- model$features
  
  
  if (!is.null(features)) {
    
    missing_features <- features[
      !features %in% colnames(new_data)
    ]
    
    
    if (length(missing_features) > 0) {
      
      stop(
        paste0(
          "The following predictor column(s) required by the model ",
          "are missing from `new_data`:\n",
          paste(
            missing_features,
            collapse = ", "
          )
        )
      )
    }
    
    
    X_new <- new_data[
      ,
      features,
      with = FALSE
    ]
    
  } else {
    
    warning(
      paste0(
        "Predictor names were not stored in this legacy model. ",
        "All numeric columns in `new_data` will be used."
      ),
      call. = FALSE
    )
    
    
    numeric_cols <- colnames(new_data)[
      vapply(
        new_data,
        is.numeric,
        logical(1)
      )
    ]
    
    
    if (length(numeric_cols) == 0) {
      
      stop(
        "No numeric predictor columns were found in `new_data`."
      )
    }
    
    
    X_new <- new_data[
      ,
      numeric_cols,
      with = FALSE
    ]
  }
  
  
  # ------------------------------------------------------------
  # Check predictor types
  # ------------------------------------------------------------
  
  non_numeric <- colnames(X_new)[
    !vapply(
      X_new,
      is.numeric,
      logical(1)
    )
  ]
  
  
  if (length(non_numeric) > 0) {
    
    stop(
      paste0(
        "The following predictor column(s) are not numeric: ",
        paste(
          non_numeric,
          collapse = ", "
        )
      )
    )
  }
  
  
  # Convert to an ordinary data.frame for model prediction.
  
  X_predict <- as.data.frame(
    X_new
  )
  
  
  # ============================================================
  # 3. Determine model type
  # ============================================================
  
  method <- model$method
  
  
  unsupervised_methods <- c(
    "one_class_svm",
    "gmm",
    "isoforest"
  )
  
  
  is_unsupervised <-
    !is.null(method) &&
    method %in% unsupervised_methods
  
  
  # ============================================================
  # 4. Apply preprocessing for unsupervised models
  # ============================================================
  
  # caret models already contain their own preprocessing.
  #
  # Unsupervised models were fitted after preprocessing outside the
  # fitted model and therefore require the stored preprocessing model
  # to be applied explicitly.
  
  if (
    is_unsupervised &&
    !is.null(model$preprocess_model)
  ) {
    
    X_predict <- stats::predict(
      model$preprocess_model,
      X_predict
    )
  }
  
  
  # ============================================================
  # 5. Generate prediction column name
  # ============================================================
  
  if (is.null(prediction_name)) {
    
    if (
      !is.null(method) &&
      method == "gmm"
    ) {
      
      prediction_name <- "cluster_predicted"
      
    } else if (!is.null(model$target)) {
      
      target_name <- tolower(
        model$target
      )
      
      
      target_name <- gsub(
        "[^A-Za-z0-9]+",
        "_",
        target_name
      )
      
      
      target_name <- gsub(
        "^_+|_+$",
        "",
        target_name
      )
      
      
      prediction_name <- paste0(
        target_name,
        "_predicted"
      )
      
    } else {
      
      prediction_name <- "prediction"
    }
  }
  
  
  # ============================================================
  # 6. Prediction
  # ============================================================
  
  start_time <- Sys.time()
  
  
  message(
    "Predicting ",
    format(
      nrow(X_predict),
      big.mark = ","
    ),
    " events..."
  )
  
  
  set.seed(
    seed
  )
  
  
  # ============================================================
  # 6a. Supervised caret models
  # ============================================================
  
  if (!is_unsupervised) {
    
    prediction_internal <- stats::predict(
      model$model,
      newdata = X_predict
    )
    
    
    prediction_internal <- as.character(
      prediction_internal
    )
    
    
    # ----------------------------------------------------------
    # Restore original biological labels
    # ----------------------------------------------------------
    
    if (!is.null(model$label_map)) {
      
      missing_labels <- unique(
        prediction_internal[
          !prediction_internal %in%
            names(model$label_map)
        ]
      )
      
      
      if (length(missing_labels) > 0) {
        
        stop(
          paste0(
            "Predicted class label(s) were not found in the model ",
            "label map: ",
            paste(
              missing_labels,
              collapse = ", "
            )
          )
        )
      }
      
      
      prediction <- unname(
        model$label_map[
          prediction_internal
        ]
      )
      
    } else {
      
      prediction <- prediction_internal
    }
    
    
    result <- data.table::data.table(
      prediction
    )
    
    
    data.table::setnames(
      result,
      "prediction",
      prediction_name
    )
    
    
    # ----------------------------------------------------------
    # Class probabilities
    # ----------------------------------------------------------
    
    if (probabilities) {
      
      set.seed(
        seed
      )
      
      
      probability_data <- tryCatch(
        
        stats::predict(
          model$model,
          newdata = X_predict,
          type = "prob"
        ),
        
        error = function(e) {
          
          stop(
            paste0(
              "Class probabilities could not be generated by this model.\n",
              "Original error: ",
              conditionMessage(e)
            )
          )
        }
      )
      
      
      probability_data <- data.table::as.data.table(
        probability_data
      )
      
      
      internal_probability_labels <- colnames(
        probability_data
      )
      
      
      if (
        !is.null(model$label_map) &&
        all(
          internal_probability_labels %in%
          names(model$label_map)
        )
      ) {
        
        probability_labels <- unname(
          model$label_map[
            internal_probability_labels
          ]
        )
        
      } else {
        
        probability_labels <-
          internal_probability_labels
      }
      
      
      probability_labels <- gsub(
        "[^A-Za-z0-9]+",
        "_",
        probability_labels
      )
      
      
      probability_names <- paste0(
        probability_labels,
        "_prob"
      )
      
      
      probability_names <- make.unique(
        probability_names,
        sep = "_"
      )
      
      
      data.table::setnames(
        probability_data,
        probability_names
      )
      
      
      result <- cbind(
        result,
        probability_data
      )
    }
  }
  
  
  # ============================================================
  # 6b. One-Class SVM
  # ============================================================
  
  if (
    !is.null(method) &&
    method == "one_class_svm"
  ) {
    
    ref_label <- model$ref_label
    
    
    other_labels <- setdiff(
      as.character(model$classes),
      ref_label
    )
    
    
    if (length(other_labels) != 1) {
      
      stop(
        "One-Class SVM prediction requires exactly two stored classes."
      )
    }
    
    
    other_label <- other_labels[1]
    
    
    prediction_raw <- stats::predict(
      model$model,
      newdata = X_predict,
      decision.values = probabilities
    )
    
    
    prediction <- ifelse(
      as.logical(prediction_raw),
      ref_label,
      other_label
    )
    
    
    result <- data.table::data.table(
      prediction
    )
    
    
    data.table::setnames(
      result,
      "prediction",
      prediction_name
    )
    
    
    if (probabilities) {
      
      decision_score <- attr(
        prediction_raw,
        "decision.values"
      )
      
      
      if (!is.null(decision_score)) {
        
        result[
          ,
          normality_score :=
            as.numeric(decision_score)
        ]
        
      } else {
        
        warning(
          paste0(
            "A decision score could not be recovered from the ",
            "One-Class SVM model."
          ),
          call. = FALSE
        )
      }
    }
  }
  
  
  # ============================================================
  # 6c. Gaussian Mixture Model
  # ============================================================
  
  if (
    !is.null(method) &&
    method == "gmm"
  ) {
    
    gmm_prediction <- stats::predict(
      model$model,
      newdata = X_predict
    )
    
    
    prediction <- as.character(
      gmm_prediction$classification
    )
    
    
    result <- data.table::data.table(
      prediction
    )
    
    
    data.table::setnames(
      result,
      "prediction",
      prediction_name
    )
    
    
    # ----------------------------------------------------------
    # Posterior cluster-membership probabilities
    # ----------------------------------------------------------
    
    if (
      probabilities &&
      !is.null(gmm_prediction$z)
    ) {
      
      probability_data <- data.table::as.data.table(
        gmm_prediction$z
      )
      
      
      probability_names <- paste0(
        "cluster_",
        seq_len(
          ncol(probability_data)
        ),
        "_prob"
      )
      
      
      data.table::setnames(
        probability_data,
        probability_names
      )
      
      
      result <- cbind(
        result,
        probability_data
      )
    }
  }
  
  
  # ============================================================
  # 6d. Isolation Forest
  # ============================================================
  
  if (
    !is.null(method) &&
    method == "isoforest"
  ) {
    
    ref_label <- model$ref_label
    
    
    other_labels <- setdiff(
      as.character(model$classes),
      ref_label
    )
    
    
    if (length(other_labels) != 1) {
      
      stop(
        "Isolation Forest prediction requires exactly two stored classes."
      )
    }
    
    
    other_label <- other_labels[1]
    
    
    # Load the Python functions distributed with AlgaeScan.
    
    isoforest_python <-
      .AlgaeScan_load_isoforest_python()
    
    
    prediction_raw <-
      isoforest_python$iso_forest_predict(
        isof_model = model$model,
        test_df = as.matrix(
          X_predict
        )
      )
    
    
    # scikit-learn IsolationForest:
    #
    #  1 = inlier
    # -1 = outlier
    
    prediction <- ifelse(
      as.integer(prediction_raw) == 1L,
      ref_label,
      other_label
    )
    
    
    result <- data.table::data.table(
      prediction
    )
    
    
    data.table::setnames(
      result,
      "prediction",
      prediction_name
    )
    
    
    # ----------------------------------------------------------
    # Normalized Isolation Forest normality score
    # ----------------------------------------------------------
    
    if (probabilities) {
      
      normality_score <-
        isoforest_python$iso_forest_predict_prob(
          isof_model = model$model,
          test_df = as.matrix(
            X_predict
          )
        )
      
      
      result[
        ,
        normality_score :=
          as.numeric(normality_score)
      ]
    }
  }
  
  
  # ============================================================
  # 7. Optionally append predictions to original data
  # ============================================================
  
  if (keep_data) {
    
    output <- data.table::copy(
      new_data
    )
    
    
    for (column in colnames(result)) {
      
      output[
        ,
        (column) := result[[column]]
      ]
    }
    
  } else {
    
    output <- result
  }
  
  
  # ============================================================
  # 8. Summary
  # ============================================================
  
  end_time <- Sys.time()
  
  
  time_taken <- difftime(
    end_time,
    start_time,
    units = "secs"
  )
  
  
  message(
    "Prediction complete in ",
    round(
      as.numeric(time_taken),
      2
    ),
    " seconds."
  )
  
  
  return(
    output
  )
}