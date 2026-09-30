#' Prepare data for AlgaeScan analysis
#'
#' Prepares imported spectral flow cytometry data for downstream
#' AlgaeScan analysis.
#'
#' The function can:
#'
#' 1. Exclude selected observations or groups.
#' 2. Attach user-defined metadata.
#' 3. Remove unwanted columns or channels.
#' 4. Split the dataset into training and testing sets.
#'
#' Metadata can contain arbitrary user-defined annotations such as
#' taxonomic class, genus, functional group, pigment group, habitat,
#' experimental condition, or other categories.
#'
#' The train/test split is performed independently within each value
#' of `stratify_by`. This reproduces the group-wise partitioning
#' strategy used in the original AlgaeScan analysis.
#'
#' @param data A `data.frame` or `data.table` containing imported
#'   spectral flow cytometry data.
#'
#' @param exclude_by Optional character string giving the column used
#'   to identify observations to exclude. Default is `NULL`.
#'
#' @param exclude_values Optional vector containing values to remove
#'   from `exclude_by`. Default is `NULL`.
#'
#' @param metadata Optional `data.frame` or `data.table` containing
#'   additional information to attach to the event data.
#'   Default is `NULL`.
#'
#' @param data_key Character string giving the identifier column in
#'   `data` used to match rows to `metadata`.
#'
#' @param metadata_key Character string giving the corresponding
#'   identifier column in `metadata`.
#'
#' @param metadata_cols Character vector containing columns from
#'   `metadata` to add to the event dataset.
#'
#' @param remove_cols Character vector containing columns to remove
#'   before partitioning. Default is `"Time"`.
#'   Use `NULL` to retain all columns.
#'
#' @param train_prop Numeric value between 0 and 1 giving the proportion
#'   of events assigned to the training dataset.
#'   Default is `0.8`.
#'
#' @param stratify_by Character string giving the column used for the
#'   stratified train/test split. Default is `"Species"`.
#'
#' @param seed Integer used for reproducible partitioning.
#'   Default is `123`.
#'
#' @param keep_data Logical. If `TRUE`, the complete prepared dataset is
#'   also returned as `data`. Keeping this additional object may require
#'   substantial memory for large cytometry datasets.
#'   Default is `FALSE`.
#'
#' @return A list containing:
#'
#' \itemize{
#'   \item `train`: training dataset.
#'   \item `test`: testing dataset.
#'   \item `data`: complete prepared dataset, only when
#'     `keep_data = TRUE`.
#' }
#'
#' @examples
#' \dontrun{
#'
#' prepared <- AlgaeScan_prepare(
#'   data = df_total_data,
#'   exclude_by = "Species",
#'   exclude_values = "RCC1502",
#'   metadata = df_species_info,
#'   data_key = "Species",
#'   metadata_key = "Roscoff Culture Collection Identifier",
#'   metadata_cols = "Class",
#'   remove_cols = "Time",
#'   train_prop = 0.8,
#'   stratify_by = "Species",
#'   seed = 123
#' )
#'
#' }
#'
#' @export
AlgaeScan_prepare <- function(
    data,
    exclude_by = NULL,
    exclude_values = NULL,
    metadata = NULL,
    data_key = NULL,
    metadata_key = NULL,
    metadata_cols = NULL,
    remove_cols = "Time",
    train_prop = 0.8,
    stratify_by = "Species",
    seed = 123,
    keep_data = FALSE
) {
  
  # ------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------
  
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
    !is.numeric(train_prop) |
    length(train_prop) != 1 |
    train_prop <= 0 |
    train_prop >= 1
  ) {
    
    stop(
      "`train_prop` must be greater than 0 and smaller than 1."
    )
  }
  
  
  # Work on one protected copy so the input supplied by
  # the user is not modified by reference.
  data <- data.table::as.data.table(data.table::copy(data))
  
  n_events_original <- nrow(data)
  
  
  # ============================================================
  # 1. Exclude selected observations
  # ============================================================
  
  if (!is.null(exclude_values)) {
    
    if (is.null(exclude_by)) {
      
      stop(
        "`exclude_by` must be specified when using `exclude_values`."
      )
    }
    
    
    if (!exclude_by %in% colnames(data)) {
      
      stop(
        paste0(
          "Column `",
          exclude_by,
          "` was not found in the dataset."
        )
      )
    }
    
    
    values_found <- exclude_values[
      exclude_values %in% data[[exclude_by]]
    ]
    
    
    values_not_found <- exclude_values[
      !exclude_values %in% data[[exclude_by]]
    ]
    
    
    if (length(values_not_found) > 0) {
      
      message(
        "Exclusion value(s) not found: ",
        paste(
          values_not_found,
          collapse = ", "
        )
      )
    }
    
    
    if (length(values_found) > 0) {
      
      n_before <- nrow(data)
      
      
      data <- data[
        !data[[exclude_by]] %in% values_found,
      ]
      
      
      n_removed <- n_before - nrow(data)
      
      
      message(
        "Excluded ",
        format(
          n_removed,
          big.mark = ","
        ),
        " events based on `",
        exclude_by,
        "`: ",
        paste(
          values_found,
          collapse = ", "
        )
      )
    }
  }
  
  
  # ============================================================
  # 2. Attach metadata
  # ============================================================
  
  if (!is.null(metadata)) {
    
    if (!is.data.frame(metadata)) {
      
      stop(
        "`metadata` must be a data.frame or data.table."
      )
    }
    
    
    if (
      is.null(data_key) |
      is.null(metadata_key) |
      is.null(metadata_cols)
    ) {
      
      stop(
        paste0(
          "When using `metadata`, specify `data_key`, ",
          "`metadata_key`, and `metadata_cols`."
        )
      )
    }
    
    
    if (!data_key %in% colnames(data)) {
      
      stop(
        paste0(
          "Column `",
          data_key,
          "` was not found in `data`."
        )
      )
    }
    
    
    if (!metadata_key %in% colnames(metadata)) {
      
      stop(
        paste0(
          "Column `",
          metadata_key,
          "` was not found in `metadata`."
        )
      )
    }
    
    
    missing_metadata_cols <- metadata_cols[
      !metadata_cols %in% colnames(metadata)
    ]
    
    
    if (length(missing_metadata_cols) > 0) {
      
      stop(
        paste0(
          "Metadata column(s) not found: ",
          paste(
            missing_metadata_cols,
            collapse = ", "
          )
        )
      )
    }
    
    
    # Metadata tables are generally small, therefore standard
    # data.frame indexing keeps this section simple and clear.
    metadata_selected <- as.data.frame(metadata)[
      ,
      c(
        metadata_key,
        metadata_cols
      ),
      drop = FALSE
    ]
    
    
    metadata_selected <- unique(
      metadata_selected
    )
    
    
    # ----------------------------------------------------------
    # Each identifier must correspond to one metadata record
    # ----------------------------------------------------------
    
    duplicated_keys <- metadata_selected[
      duplicated(
        metadata_selected[[metadata_key]]
      ) |
        duplicated(
          metadata_selected[[metadata_key]],
          fromLast = TRUE
        ),
      ,
      drop = FALSE
    ]
    
    
    if (nrow(duplicated_keys) > 0) {
      
      problem_values <- unique(
        duplicated_keys[[metadata_key]]
      )
      
      
      stop(
        paste0(
          "Some values in `",
          metadata_key,
          "` map to more than one metadata record: ",
          paste(
            problem_values,
            collapse = ", "
          )
        )
      )
    }
    
    
    # ----------------------------------------------------------
    # Match metadata without changing event order
    # ----------------------------------------------------------
    
    metadata_index <- match(
      data[[data_key]],
      metadata_selected[[metadata_key]]
    )
    
    
    for (column in metadata_cols) {
      
      data[[column]] <-
        metadata_selected[[column]][metadata_index]
    }
    
    
    unmatched_values <- unique(
      data[[data_key]][is.na(metadata_index)]
    )
    
    
    if (length(unmatched_values) > 0) {
      
      message(
        "No metadata found for ",
        length(unmatched_values),
        " value(s) in `",
        data_key,
        "`. Metadata set to NA."
      )
    }
    
    
    message(
      "Added metadata column(s): ",
      paste(
        metadata_cols,
        collapse = ", "
      )
    )
  }
  
  
  # ============================================================
  # 3. Remove selected columns
  # ============================================================
  
  if (!is.null(remove_cols)) {
    
    cols_found <- remove_cols[
      remove_cols %in% colnames(data)
    ]
    
    
    cols_not_found <- remove_cols[
      !remove_cols %in% colnames(data)
    ]
    
    
    if (length(cols_not_found) > 0) {
      
      message(
        "Column(s) not found and therefore not removed: ",
        paste(
          cols_not_found,
          collapse = ", "
        )
      )
    }
    
    
    if (length(cols_found) > 0) {
      
      data[
        ,
        (cols_found) := NULL
      ]
      
      
      message(
        "Removed column(s): ",
        paste(
          cols_found,
          collapse = ", "
        )
      )
    }
  }
  
  
  # ============================================================
  # 4. Check stratification variable
  # ============================================================
  
  if (!stratify_by %in% colnames(data)) {
    
    stop(
      paste0(
        "Column `",
        stratify_by,
        "` was not found in the prepared dataset."
      )
    )
  }
  
  
  if (any(is.na(data[[stratify_by]]))) {
    
    stop(
      paste0(
        "`",
        stratify_by,
        "` contains missing values."
      )
    )
  }
  
  
  # ============================================================
  # 5. Stratified train/test split
  #
  # Memory-efficient version:
  # only row indices are stored while looping through groups.
  #
  # The seed is reset inside each group exactly as in the
  # original partition_data() implementation.
  # ============================================================
  
  groups <- unique(
    data[[stratify_by]]
  )
  
  
  train_indices <- vector(
    "list",
    length(groups)
  )
  
  
  test_indices <- vector(
    "list",
    length(groups)
  )
  
  
  names(train_indices) <- groups
  names(test_indices) <- groups
  
  
  for (i in seq_along(groups)) {
    
    group <- groups[i]
    
    
    # Global row positions belonging to this group
    group_rows <- which(
      data[[stratify_by]] == group
    )
    
    
    set.seed(seed)
    
    
    # createDataPartition returns row positions relative
    # to the current group
    train_local <- caret::createDataPartition(
      y = data[[stratify_by]][group_rows],
      p = train_prop,
      list = TRUE
    )[[1]]
    
    
    train_indices[[i]] <-
      group_rows[train_local]
    
    
    test_indices[[i]] <-
      group_rows[-train_local]
  }
  
  
  # Combine indices
  train_indices <- unlist(
    train_indices,
    use.names = FALSE
  )
  
  
  test_indices <- unlist(
    test_indices,
    use.names = FALSE
  )
  
  
  # Only two full datasets are now created
  train_data <- data[
    train_indices,
  ]
  
  
  test_data <- data[
    test_indices,
  ]
  
  
  # ============================================================
  # Summary
  # ============================================================
  
  message("")
  message("AlgaeScan preparation complete.")
  
  message(
    "Events before preparation: ",
    format(
      n_events_original,
      big.mark = ","
    )
  )
  
  message(
    "Prepared events: ",
    format(
      nrow(data),
      big.mark = ","
    )
  )
  
  message(
    "Training events: ",
    format(
      nrow(train_data),
      big.mark = ","
    )
  )
  
  message(
    "Testing events: ",
    format(
      nrow(test_data),
      big.mark = ","
    )
  )
  
  message(
    "Split: ",
    train_prop * 100,
    "/",
    (1 - train_prop) * 100,
    " by `",
    stratify_by,
    "`"
  )
  
  
  # ============================================================
  # Return
  # ============================================================
  
  if (keep_data) {
    
    return(
      list(
        data = data,
        train = train_data,
        test = test_data
      )
    )
    
  } else {
    
    # Do not retain the third full dataset in the returned object
    return(
      list(
        train = train_data,
        test = test_data
      )
    )
  }
}


#' Downsample events by group
#'
#' Randomly downsamples events within groups using either a fixed number
#' of events per group or a proportion of events per group.
#'
#' The same downsampling value can be applied to every group by supplying
#' a single value. Different values can be applied to different groups
#' using a named vector.
#'
#' The relative order of selected rows in the original dataset is
#' preserved.
#'
#' @param data A data.frame or data.table containing the data to
#'   downsample.
#'
#' @param n Number of events to retain per group. Can be either:
#'   \itemize{
#'     \item a single number applied to every group;
#'     \item a named numeric vector giving a different number for each
#'     group.
#'   }
#'   Default is `300`. Do not use together with `prop`.
#'
#' @param by Character string giving the column used to define groups.
#'   For example `"Species"` or `"type_predicted"`.
#'
#' @param seed Random seed used for sampling. Default is `123`.
#'
#' @param prop Proportion of events to retain per group. Can be either:
#'   \itemize{
#'     \item a single value between 0 and 1 applied to every group;
#'     \item a named numeric vector giving a different proportion for
#'     each group.
#'   }
#'   Do not use together with `n`.
#'
#' @return A data.table containing the selected events.
#'
#' @examples
#' \dontrun{
#'
#' # Keep up to 300 events per strain
#' AlgaeScan_downsample(
#'   data = data,
#'   n = 300,
#'   by = "Species",
#'   seed = 123
#' )
#'
#' # Keep 10% of events from every strain
#' AlgaeScan_downsample(
#'   data = data,
#'   prop = 0.10,
#'   by = "Species",
#'   seed = 123
#' )
#'
#' # Use different proportions for different predicted classes
#' AlgaeScan_downsample(
#'   data = data,
#'   prop = c(
#'     "algae" = 0.001,
#'     "non-algae" = 0.10
#'   ),
#'   by = "type_predicted",
#'   seed = 123
#' )
#'
#' }
#'
#' @export
AlgaeScan_downsample <- function(
    data,
    n = 300,
    by,
    seed = 123,
    prop = NULL
) {
  
  # ============================================================
  # 1. Check input data
  # ============================================================
  
  if (!is.data.frame(data)) {
    
    stop(
      "`data` must be a data.frame or data.table."
    )
  }
  
  
  if (nrow(data) == 0) {
    
    stop(
      "`data` contains no rows."
    )
  }
  
  
  data <- data.table::as.data.table(
    data.table::copy(data)
  )
  
  
  # ============================================================
  # 2. Check grouping variable
  # ============================================================
  
  if (
    !is.character(by) ||
    length(by) != 1
  ) {
    
    stop(
      "`by` must be a single character string."
    )
  }
  
  
  if (!by %in% colnames(data)) {
    
    stop(
      paste0(
        "Grouping column `",
        by,
        "` was not found."
      )
    )
  }
  
  
  if (anyNA(data[[by]])) {
    
    stop(
      paste0(
        "Grouping column `",
        by,
        "` contains missing values."
      )
    )
  }
  
  
  # ============================================================
  # 3. Determine downsampling mode
  # ============================================================
  
  # If prop is supplied and n was not explicitly supplied,
  # use proportional downsampling.
  #
  # This keeps the existing default n = 300 behavior while
  # allowing:
  #
  # AlgaeScan_downsample(
  #   data = data,
  #   prop = 0.10,
  #   by = "Species"
  # )
  
  if (!is.null(prop)) {
    
    if (
      !missing(n) &&
      !is.null(n)
    ) {
      
      stop(
        "Use either `n` or `prop`, not both."
      )
    }
    
    
    mode <- "prop"
    
  } else {
    
    if (is.null(n)) {
      
      stop(
        "Either `n` or `prop` must be supplied."
      )
    }
    
    
    mode <- "n"
  }
  
  
  # ============================================================
  # 4. Split row indices by group
  # ============================================================
  
  groups <- split(
    seq_len(nrow(data)),
    as.character(data[[by]])
  )
  
  
  group_names <- names(
    groups
  )
  
  
  # ============================================================
  # 5. Prepare number-based downsampling
  # ============================================================
  
  if (mode == "n") {
    
    if (
      !is.numeric(n) ||
      anyNA(n) ||
      any(!is.finite(n))
    ) {
      
      stop(
        "`n` must contain numeric, finite values."
      )
    }
    
    
    if (any(n < 1)) {
      
      stop(
        "`n` values must be at least 1."
      )
    }
    
    
    if (
      any(
        abs(n - round(n)) >
        sqrt(.Machine$double.eps)
      )
    ) {
      
      stop(
        "`n` must contain whole numbers."
      )
    }
    
    
    # ----------------------------------------------------------
    # Same n for every group
    # ----------------------------------------------------------
    
    if (length(n) == 1) {
      
      values_by_group <- rep(
        as.integer(n),
        length(group_names)
      )
      
      
      names(
        values_by_group
      ) <- group_names
      
      
    } else {
      
      # --------------------------------------------------------
      # Different n for different groups
      # --------------------------------------------------------
      
      if (
        is.null(names(n)) ||
        any(names(n) == "")
      ) {
        
        stop(
          paste0(
            "When multiple values are supplied to `n`, ",
            "`n` must be a named vector."
          )
        )
      }
      
      
      if (anyDuplicated(names(n))) {
        
        stop(
          "`n` contains duplicated group names."
        )
      }
      
      
      missing_groups <- setdiff(
        group_names,
        names(n)
      )
      
      
      if (length(missing_groups) > 0) {
        
        stop(
          paste0(
            "No `n` value was supplied for: ",
            paste(
              missing_groups,
              collapse = ", "
            )
          )
        )
      }
      
      
      values_by_group <- as.integer(
        n[
          group_names
        ]
      )
      
      
      names(
        values_by_group
      ) <- group_names
    }
  }
  
  
  # ============================================================
  # 6. Prepare proportion-based downsampling
  # ============================================================
  
  if (mode == "prop") {
    
    if (
      !is.numeric(prop) ||
      anyNA(prop) ||
      any(!is.finite(prop))
    ) {
      
      stop(
        "`prop` must contain numeric, finite values."
      )
    }
    
    
    if (
      any(prop <= 0) ||
      any(prop > 1)
    ) {
      
      stop(
        "`prop` values must be greater than 0 and less than or equal to 1."
      )
    }
    
    
    # ----------------------------------------------------------
    # Same proportion for every group
    # ----------------------------------------------------------
    
    if (length(prop) == 1) {
      
      values_by_group <- rep(
        prop,
        length(group_names)
      )
      
      
      names(
        values_by_group
      ) <- group_names
      
      
    } else {
      
      # --------------------------------------------------------
      # Different proportions for different groups
      # --------------------------------------------------------
      
      if (
        is.null(names(prop)) ||
        any(names(prop) == "")
      ) {
        
        stop(
          paste0(
            "When multiple values are supplied to `prop`, ",
            "`prop` must be a named vector."
          )
        )
      }
      
      
      if (anyDuplicated(names(prop))) {
        
        stop(
          "`prop` contains duplicated group names."
        )
      }
      
      
      missing_groups <- setdiff(
        group_names,
        names(prop)
      )
      
      
      if (length(missing_groups) > 0) {
        
        stop(
          paste0(
            "No `prop` value was supplied for: ",
            paste(
              missing_groups,
              collapse = ", "
            )
          )
        )
      }
      
      
      values_by_group <- prop[
        group_names
      ]
    }
  }
  
  
  # ============================================================
  # 7. Downsample each group
  # ============================================================
  
  set.seed(seed)
  
  
  selected_rows <- vector(
    mode = "list",
    length = length(groups)
  )
  
  
  for (i in seq_along(groups)) {
    
    rows_i <- groups[[i]]
    
    
    # ----------------------------------------------------------
    # Number of events
    # ----------------------------------------------------------
    
    if (mode == "n") {
      
      n_take <- min(
        values_by_group[i],
        length(rows_i)
      )
    }
    
    
    # ----------------------------------------------------------
    # Proportion of events
    # ----------------------------------------------------------
    
    if (mode == "prop") {
      
      n_take <- ceiling(
        length(rows_i) *
          values_by_group[i]
      )
      
      
      n_take <- min(
        n_take,
        length(rows_i)
      )
    }
    
    
    # ----------------------------------------------------------
    # Sample rows
    # ----------------------------------------------------------
    
    if (n_take == length(rows_i)) {
      
      selected_rows[[i]] <- rows_i
      
    } else {
      
      selected_rows[[i]] <- sample(
        rows_i,
        size = n_take,
        replace = FALSE
      )
    }
  }
  
  
  # ============================================================
  # 8. Restore original row order
  # ============================================================
  
  selected_rows <- sort(
    unlist(
      selected_rows,
      use.names = FALSE
    )
  )
  
  
  output <- data[
    selected_rows
  ]
  
  
  # ============================================================
  # 9. Return
  # ============================================================
  
  return(
    output
  )
}