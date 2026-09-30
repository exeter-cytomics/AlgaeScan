#' Perform UMAP dimensionality reduction
#'
#' Performs UMAP dimensionality reduction on spectral flow cytometry data.
#' Numeric columns are automatically selected as candidate features.
#' Columns can be removed using regular-expression patterns.
#'
#' Optional preprocessing includes area-channel selection, asinh
#' transformation, unit vector normalization, and PCA.
#'
#' An existing UMAP and PCA model can also be supplied to project new
#' data into an existing UMAP space.
#'
#' @param data A data.frame or data.table containing spectral flow
#'   cytometry data.
#'
#' @param remove_cols Optional character vector containing regular
#'   expressions identifying columns to remove before UMAP.
#'   Default is `c("^Time$", "_prob$")`.
#'
#' @param area_only Logical. If `TRUE`, retain only channels ending
#'   in `"-A"`. Default is `FALSE`.
#'
#' @param n_neighbors Number of nearest neighbours used by UMAP.
#'   Default is `5`.
#'
#' @param min_dist Minimum distance used by UMAP.
#'   Default is `0.05`.
#'
#' @param n_epochs Number of UMAP training epochs.
#'   Default is `500`.
#'
#' @param spread Spread parameter used by UMAP.
#'   Default is `1`.
#'
#' @param metric Distance metric used by UMAP.
#'   Default is `"euclidean"`.
#'
#' @param asinh_transform Logical. If `TRUE`, apply an asinh
#'   transformation before UMAP. Default is `FALSE`.
#'
#' @param cofactor_asinh Cofactor used for the asinh transformation.
#'   Default is `150`.
#'
#' @param apply_pca Logical. If `TRUE`, perform PCA before UMAP.
#'   Default is `FALSE`.
#'
#' @param pca_components Number of principal components retained before
#'   UMAP. Default is `30`.
#'
#' @param init Initialization method used by UMAP.
#'   Default is `"spectral"`.
#'
#' @param uvn Logical. If `TRUE`, perform unit vector normalization.
#'   Default is `FALSE`.
#'
#' @param umap_model Optional previously fitted UMAP model. When supplied,
#'   the new data are projected into the existing UMAP space.
#'
#' @param pca_model Optional previously fitted PCA model. Required when
#'   projecting data with `apply_pca = TRUE`.
#'
#' @param seed Random seed.
#'   Default is `123`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `umap_df`: UMAP coordinates and available metadata.
#'     \item `umap_model`: fitted UMAP model.
#'     \item `pca_model`: fitted PCA model, when PCA is used.
#'     \item `features`: predictor columns used for dimensionality reduction.
#'     \item `settings`: UMAP and preprocessing settings.
#'   }
#'
#' @examples
#' \dontrun{
#'
#' # Remove FSC and SSC channels and use fluorescence Area channels
#' out <- AlgaeScan_umap(
#'   data = data,
#'   remove_cols = c("^Time$", "_prob$", "FSC", "SSC"),
#'   area_only = TRUE
#' )
#'
#' }
#'
#' @export
AlgaeScan_umap <- function(
    data,
    remove_cols = c("^Time$", "_prob$"),
    area_only = FALSE,
    n_neighbors = 5,
    min_dist = 0.05,
    n_epochs = 500,
    spread = 1,
    metric = "euclidean",
    asinh_transform = FALSE,
    cofactor_asinh = 150,
    apply_pca = FALSE,
    pca_components = 30,
    init = "spectral",
    uvn = FALSE,
    umap_model = NULL,
    pca_model = NULL,
    seed = 123
) {
  
  # ============================================================
  # 1. Check input
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
  
  if (!is.null(remove_cols) && !is.character(remove_cols)) {
    stop(
      "`remove_cols` must be NULL or a character vector."
    )
  }
  
  if (!is.numeric(seed) || length(seed) != 1) {
    stop(
      "`seed` must be a single numeric value."
    )
  }
  
  if (
    !is.numeric(pca_components) ||
    length(pca_components) != 1 ||
    pca_components < 1
  ) {
    stop(
      "`pca_components` must be a single positive number."
    )
  }
  
  data <- data.table::as.data.table(
    data
  )
  
  
  # ============================================================
  # 2. Select numeric predictor columns
  # ============================================================
  
  message("Checking predictor columns...")
  
  numeric_cols <- colnames(data)[
    sapply(
      data,
      is.numeric
    )
  ]
  
  if (length(numeric_cols) == 0) {
    stop(
      "No numeric columns were found in `data`."
    )
  }
  
  X <- data[
    ,
    numeric_cols,
    with = FALSE
  ]
  
  
  # ============================================================
  # 3. Remove unwanted columns
  # ============================================================
  
  if (!is.null(remove_cols)) {
    
    cols_to_remove <- unique(
      unlist(
        lapply(
          remove_cols,
          function(pattern) {
            
            # Exact column-name match
            if (pattern %in% colnames(X)) {
              
              return(pattern)
              
            } else {
              
              # Otherwise interpret as regex
              return(
                grep(
                  pattern,
                  colnames(X),
                  value = TRUE,
                  perl = TRUE
                )
              )
            }
          }
        )
      )
    )
    
    
    if (length(cols_to_remove) > 0) {
      
      message(
        "Removing ",
        length(cols_to_remove),
        " column(s): ",
        paste(
          cols_to_remove,
          collapse = ", "
        )
      )
      
      
      X[
        ,
        (cols_to_remove) := NULL
      ]
    }
  }
  
  # ============================================================
  # 4. Optionally retain Area channels only
  # ============================================================
  
  if (area_only) {
    
    message(
      "Keeping Area channels..."
    )
    
    area_cols <- grep(
      "-A$",
      colnames(X),
      value = TRUE
    )
    
    
    if (length(area_cols) == 0) {
      
      stop(
        "No channels ending in `-A` were found."
      )
    }
    
    
    X <- X[
      ,
      area_cols,
      with = FALSE
    ]
  }
  
  
  # ============================================================
  # 5. Check selected predictors
  # ============================================================
  
  if (ncol(X) == 0) {
    
    stop(
      "No predictor columns remain after column filtering."
    )
  }
  
  
  non_numeric <- colnames(X)[
    !sapply(
      X,
      is.numeric
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
  
  
  feature_names <- colnames(
    X
  )
  
  
  message(
    "Using ",
    length(feature_names),
    " features for UMAP."
  )
  
  
  X <- as.matrix(
    X
  )
  
  
  # ============================================================
  # 6. Asinh transformation
  # ============================================================
  
  if (asinh_transform) {
    
    message(
      "Applying asinh transformation..."
    )
    
    X <- asinh(
      X / cofactor_asinh
    )
  }
  
  
  # ============================================================
  # 7. Unit vector normalization
  # ============================================================
  
  if (uvn) {
    
    message(
      "Applying unit vector normalization..."
    )
    
    row_norms <- sqrt(
      rowSums(
        X^2
      )
    )
    
    
    # Prevent division by zero
    row_norms[
      row_norms == 0
    ] <- 1
    
    
    X <- X / row_norms
  }
  
  
  # ============================================================
  # 8. PCA + UMAP
  # ============================================================
  
  set.seed(seed)
  
  
  if (is.null(umap_model)) {
    
    # ==========================================================
    # Fit new embedding
    # ==========================================================
    
    if (apply_pca) {
      
      message(
        "Applying PCA..."
      )
      
      
      max_components <- min(
        ncol(X),
        nrow(X) - 1
      )
      
      
      if (pca_components > max_components) {
        
        stop(
          paste0(
            "`pca_components` cannot exceed ",
            max_components,
            " for this dataset."
          )
        )
      }
      
      
      pca_model <- stats::prcomp(
        X,
        rank. = pca_components,
        center = TRUE,
        scale. = TRUE
      )
      
      
      X_umap <- pca_model$x
      
    } else {
      
      X_umap <- X
    }
    
    
    # ----------------------------------------------------------
    # Fit UMAP
    # ----------------------------------------------------------
    
    message(
      "Applying UMAP..."
    )
    
    
    set.seed(seed)
    
    
    umap_result <- umap::umap(
      X_umap,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      metric = metric,
      init = init,
      n_epochs = n_epochs,
      spread = spread,
      random_state = seed
    )
    
    
    embedding <- umap_result$layout
    
    fitted_umap_model <- umap_result
    
    
  } else {
    
    # ==========================================================
    # Project data into existing UMAP space
    # ==========================================================
    
    if (apply_pca) {
      
      if (is.null(pca_model)) {
        
        stop(
          paste0(
            "`pca_model` must be supplied when using an existing ",
            "UMAP model with `apply_pca = TRUE`."
          )
        )
      }
      
      
      # --------------------------------------------------------
      # Recover exact features and order from PCA model
      # --------------------------------------------------------
      
      pca_features <- rownames(
        pca_model$rotation
      )
      
      
      missing_features <- pca_features[
        !pca_features %in% feature_names
      ]
      
      
      if (length(missing_features) > 0) {
        
        stop(
          paste0(
            "The following predictor column(s) required by the PCA model ",
            "are missing:\n",
            paste(
              missing_features,
              collapse = ", "
            )
          )
        )
      }
      
      
      feature_positions <- match(
        pca_features,
        feature_names
      )
      
      
      X <- X[
        ,
        feature_positions,
        drop = FALSE
      ]
      
      
      message(
        "Applying existing PCA model..."
      )
      
      
      X_umap <- stats::predict(
        pca_model,
        X
      )
      
      
    } else {
      
      X_umap <- X
    }
    
    
    # ----------------------------------------------------------
    # Apply existing UMAP
    # ----------------------------------------------------------
    
    message(
      "Applying existing UMAP model..."
    )
    
    
    set.seed(seed)
    
    
    embedding <- predict(
      umap_model,
      X_umap
    )
    
    
    fitted_umap_model <- umap_model
  }
  
  
  # ============================================================
  # 9. Build output table
  # ============================================================
  
  message(
    "Extracting UMAP results..."
  )
  
  
  umap_df <- data.table::data.table(
    UMAP1 = embedding[, 1],
    UMAP2 = embedding[, 2]
  )
  
  
  # ============================================================
  # 10. Add metadata
  # ============================================================
  
  metadata_cols <- c(
    "Species",
    "measure",
    "type",
    "type_predicted",
    "type_predicted_prob_thr",
    "species_predicted",
    "Species_with_noise",
    "species_predicted_with_noise",
    "Class",
    "class"
  )
  
  
  metadata_cols <- metadata_cols[
    metadata_cols %in% colnames(data)
  ]
  
  
  # ------------------------------------------------------------
  # Keep probability columns as metadata
  # ------------------------------------------------------------
  
  probability_cols <- grep(
    "_prob$",
    colnames(data),
    value = TRUE
  )
  
  
  # ------------------------------------------------------------
  # Keep batch information
  # ------------------------------------------------------------
  
  batch_cols <- grep(
    "batch",
    colnames(data),
    value = TRUE
  )
  
  
  metadata_cols <- unique(
    c(
      metadata_cols,
      probability_cols,
      batch_cols
    )
  )
  
  
  if (length(metadata_cols) > 0) {
    
    metadata <- data[
      ,
      metadata_cols,
      with = FALSE
    ]
    
    
    if (
      "Class" %in% colnames(metadata) &&
      !"class" %in% colnames(metadata)
    ) {
      
      data.table::setnames(
        metadata,
        "Class",
        "class"
      )
    }
    
    
    umap_df <- cbind(
      umap_df,
      metadata
    )
  }
  
  
  # ============================================================
  # 11. Return results
  # ============================================================
  
  return(
    list(
      
      umap_df = umap_df,
      
      umap_model = fitted_umap_model,
      
      pca_model = pca_model,
      
      features = feature_names,
      
      settings = list(
        
        remove_cols = remove_cols,
        
        area_only = area_only,
        
        n_neighbors = n_neighbors,
        
        min_dist = min_dist,
        
        n_epochs = n_epochs,
        
        spread = spread,
        
        metric = metric,
        
        asinh_transform = asinh_transform,
        
        cofactor_asinh = cofactor_asinh,
        
        apply_pca = apply_pca,
        
        pca_components = pca_components,
        
        init = init,
        
        uvn = uvn,
        
        seed = seed
      )
    )
  )
}