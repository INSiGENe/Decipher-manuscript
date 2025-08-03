#' Add a Name Column to a Data Frame
#'
#' This function adds a new column to an existing data frame where all entries in
#' the new column are set to the specified name. This can be useful for tracking the
#' origin of data after merging several data frames.
#'
#' @param df A data frame to which the new column will be added.
#' @param name A character string that will be used to fill the new column.
#'
#' @return The modified data frame with an additional column named 'DataFrameName'
#' which contains the same value specified by the `name` parameter for all rows.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'your_dataframe' is your existing data frame and you want to
#'   # label it as "sample_data"
#'   labeled_df <- add_name_column(your_dataframe, "sample_data")
#'   print(labeled_df)
#' }
#'
#' @export
add_name_column <- function(df, name) {
  df$DataFrameName <- name
  return(df)
}

#' Calculate Absolute Maximum Value
#'
#' Computes the maximum value in absolute terms from a numeric vector. This
#' function is useful for identifying the largest magnitude value, regardless
#' of its sign, within a given set of numbers.
#'
#' @param row A numeric vector for which the absolute maximum value is sought.
#'
#' @return A single numeric value representing the maximum absolute value
#' found in the input vector.
#'
#' @examples
#' \dontrun{
#'   sample_vector <- c(-10, 2, 3, -4, 5)
#'   max_val <- get_abs_max(sample_vector)
#'   print(max_val) # Outputs: 10
#' }
#'
#' @export
get_abs_max <- function(row) {
  max(abs(row))
}

#' Filter Rows with Class "Real"
#'
#' This function filters rows in a dataframe where the `class` column is equal to "real".
#' It's particularly useful for datasets where observations are categorized by class,
#' allowing for the isolation of real instances from others.
#'
#' @param df A dataframe containing at least one column named `class`, with mixed class types.
#'
#' @return A dataframe consisting only of rows where the `class` column has the value "real".
#'
#' @examples
#' \dontrun{
#'   # Assuming 'data_df' is your dataframe that includes a 'class' column
#'   real_df <- filter_real(data_df)
#'   print(real_df)
#' }
#'
#' @importFrom dplyr filter
#' @export
filter_real <- function(df) {
  df %>% filter(class == "real")
}


#' Convert Specific Text Patterns in Strings
#'
#' This function takes a vector of strings and replaces specific text patterns
#' according to the following rules:
#' - Replaces "_plus_" with "+"
#' - Replaces "_minus_" with "-"
#' - Replaces "monocytes" with " Mono"
#'
#' @param text_vec A character vector containing the text to be processed.
#'
#' @return A character vector with the specified text patterns replaced.
#'
#' @examples
#' text_vec <- c("T_cells_plus_", "B_cells_minus_", "monocytes")
#' converted_text <- convert_text_patterns(text_vec)
#' print(converted_text)
#'
#' @export
convert_text_patterns <- function(text_vec) {
  text_vec <- gsub("_plus_", "+", text_vec)
  text_vec <- gsub("_minus_", "-", text_vec)
  text_vec <- gsub("monocytes", "Mono", text_vec)
  text_vec <- gsub("CD14\\+BDCA1\\+PD-L1\\+cells", "C8", text_vec)

  return(text_vec)
}


#' Check and Adapt Decipher Dataframe Format
#'
#' Examines the first dataframe within a list of dataframes to determine if it uses an
#' older version of Decipher naming conventions (specifically, checks for the presence
#' of "delta.pagoda" column name). If this older convention is detected, the function
#' updates the column name to the newer version ("deltaPagoda") across all dataframes in
#' the list. It prints a message indicating whether the data was updated or is already
#' in the expected format.
#'
#' @param regulon_deltas_by_cluster A list of dataframes, each representing regulon
#' delta values by cluster. Expected to possibly contain a column named "delta.pagoda",
#' which will be renamed to "deltaPagoda".
#'
#' @return The same list of dataframes passed as input, potentially with column names
#' updated from "delta.pagoda" to "deltaPagoda".
#'
#' @examples
#' \dontrun{
#'   # Assuming 'regulon_deltas' is your list of dataframes with regulon delta values
#'   updated_deltas <- checkAndAdaptDecipherVersion(regulon_deltas)
#' }
#'
#' @importFrom dplyr rename
#' @export
checkAndAdaptDecipherVersion <- function(regulon_deltas_by_cluster){
  if(colnames(regulon_deltas_by_cluster[[1]])[1] %in% c("delta.pagoda")){
    print("this appears to be Decipher before refactoring")
    regulon_deltas_by_cluster <- lapply(regulon_deltas_by_cluster, function(df) {
      df %>% rename(deltaPagoda = delta.pagoda)
    })}
  else{print("this Decipher version is up to date")}
  return(regulon_deltas_by_cluster)
}

#' Replace NA values with 0 in a vector
#'
#' This function takes a vector and replaces all NA values with 0.
#' The modified vector with 0s replacing the NAs is then returned.
#'
#' @param vector A numeric vector that may contain NA values.
#' @return A numeric vector where all NA values have been replaced with 0.
#' @examples
#' sample_vector <- c(1, NA, 3, NA, 5)
#' replaceNAw0(sample_vector)
#' @export
replaceNAw0 <- function(vector){
  vector[is.na(vector)] <- 0
  return(vector)
}


#' Append a new entry to a data frame
#'
#' This function takes an existing data frame and a new entry (as a row) and appends the new entry to the original data frame using `rbind`.
#' The resulting data frame, with the new entry added, is then returned.
#'
#' @param original_df A data frame to which the new entry will be added.
#' @param new_entry A list or a data frame row that represents the new entry to be added.
#' @return A data frame with the new entry appended to the original data frame.
#' @examples
#' original_df <- data.frame(Name = c("Alice", "Bob"), Age = c(25, 30))
#' new_entry <- data.frame(Name = "Charlie", Age = 35)
#' addEntryToDF(original_df, new_entry)
#' @export
addEntryToDF <- function(original_df,new_entry){
  return(rbind(original_df,new_entry))
}

#' Clean special characters in string names
#'
#' This function modifies strings by replacing specific patterns with more conventional characters.
#' It replaces '_minus_' with '-' and '_plus_' with '+'. This can be particularly useful in cleaning up variable names or labels that have been encoded with special characters to adhere to programming language syntax constraints.
#'
#' @param this_string A single string or character vector where patterns need to be replaced.
#' @return A string or character vector with the specified patterns replaced by conventional characters.
#' @examples
#' clean_names("rate_minus_inflation_plus_growth")
#' @importFrom stringr str_replace_all
#' @export
clean_names <- function(this_string){
  this_string <- stringr::str_replace_all(this_string, pattern="_minus_", replacement="-")
  this_string <- stringr::str_replace_all(this_string, pattern="_plus_", replacement="+")
  return(this_string)
}

#' Filter Data Based on Threshold
#'
#' This function filters a dataset based on two thresholds. It selects rows that meet the main threshold, and if none do, it selects rows with the fallback threshold.
#'
#' @param data A data frame containing the dataset to be filtered.
#' @param threshold_main The main threshold value used for filtering.
#' @param threshold_fallback The fallback threshold value used if no rows meet the main threshold.
#'
#' @return A subset of the input data frame containing rows that meet either the main threshold or the fallback threshold.
#'
#' @examples
#' my_data <- data.frame(
#'   ID = c(1, 2, 3, 4, 5),
#'   Value = c(10, 20, 30, 40, 50),
#'   Threshold = c(15, 25, 35, 45, 55)
#' )
#' filter_threshold(my_data, 30, 25)
#'
#' @export
filter_threshold <- function(data, threshold_main, threshold_fallback) {
  # If any row meets the main threshold, select those
  if (any(data$threshold == threshold_main)) {
    return(data[data$threshold == threshold_main, ])
  } else {
    # Otherwise, select rows with the fallback threshold
    return(data[data$threshold == threshold_fallback, ])
  }
}

#' Assign Family Colors to Each Regulon
#'
#' This function takes a matrix of the top 30 regulons and their delta values, and assigns a specific color
#' to each transcription factor (TF) based on its family. The five TF families included are NFKB, STAT, IRF,
#' AP1, and NFAT. Transcription factors that are not in these families are colored gray.
#'
#' @param top_30_regulons_delta_matrix A matrix with regulons as rows. The names of the regulons must be set
#'        as the row names of the matrix. It is assumed that the matrix contains the top 30 regulons based on
#'        delta values.
#'
#' @return A named vector of colors corresponding to the transcription factors in the input matrix. The name of
#'         each element in the vector corresponds to a regulon, and the value is the color assigned based on
#'         the transcription factor family.
#'
#' @examples
#' # Assuming delta_matrix is a matrix with row names that include transcription factors
#' tf_colors <- assign_tf_family_colors_to_each_regulon(delta_matrix)
#' print(tf_colors)
#'
#' @export
assign_tf_family_colors_to_each_regulon <- function(top_30_regulons_delta_matrix){
  families <- list(
    NFKB = c("NFKB1","NFKB2","REL","RELA","RELB"),
    STAT = c("STAT1","STAT2","STAT3","STAT4","STAT5B"),
    IRF = c("IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8","IRF9"),
    AP1 = c("JUN","JUND","JUNB","FOS","FOSL2","ATF4","ATF5","ATF3","BATF","BATF3","ATF1","ATF2","MAFK","MAF","MAFG","MAFB","MAFF"),
    NFAT = c("NFATC1","NFAT5","NFATC2")
  )

  # Create a named vector for TF colors
  tf_colors <- rep("gray", length(rownames(top_30_regulons_delta_matrix)))
  names(tf_colors) <- rownames(top_30_regulons_delta_matrix)

  # Assign colors to TFs based on their family
  tf_colors[names(tf_colors) %in% unlist(families$IRF)] <- "red"
  tf_colors[names(tf_colors) %in% unlist(families$AP1)] <- "blue"
  tf_colors[names(tf_colors) %in% unlist(families$NFKB)] <- "green"
  tf_colors[names(tf_colors) %in% unlist(families$STAT)] <- "purple"
  tf_colors[names(tf_colors) %in% unlist(families$NFAT)] <- "yellow"

  return(tf_colors)
}


#' Safely load a CSV or RDS file
#'
#' Attempts to read a file based on its extension.  
#' If the file does not exist or an error occurs during loading,  
#' a message is emitted and `NULL` is returned.
#'
#' @param filepath Character. Path to the file to load.
#' @param ... Additional arguments passed to `read.csv()` when loading CSV files.
#'
#' @return The contents of the file (usually a data frame or R object), or `NULL` if the file  
#'   is missing, unsupported, or an error is thrown.
#'
#'
#' @export
safe_load <- function(filepath, ...) {
  if (!file.exists(filepath)) {
    message("Skipping missing file: ", filepath)
    return(NULL)
  }

  ext <- tools::file_ext(filepath)

  tryCatch({
    switch(ext,
      "rds" = readRDS(filepath),
      "csv" = read.csv(filepath, ...),
      stop(paste("Unsupported file type:", ext))
    )
  }, error = function(e) {
    message("Error loading file: ", filepath, " — ", e$message)
    return(NULL)
  })
}



#' Conditionally add an element to a list
#'
#' Adds \`value\` into the provided list under \`name\` if \`value\` is not \`NULL\`.  
#' If \`value\` is a data frame, only the columns \`sender\`, \`receiver\`,  
#' \`interaction\`, \`prioritization_score\`, and \`scaled_score\` are retained.
#'
#' @param lst A named list to which the new element may be added.
#' @param name Character. The name under which to store \`value\` in \`lst\`.
#' @param value The object to add. If it is \`NULL\`, nothing is added.  
#'   If it is a data frame, only the specified columns are kept.
#'
#' @return The updated list, with \`value\` added under \`name\` if it was non-NULL.
#'
#' @importFrom dplyr select
#' @export
add_if_not_null <- function(lst, name, value) {
  if (!is.null(value)) {
    if (is.data.frame(value)) {
      lst[[name]] <- value %>% select(sender, receiver, interaction, prioritization_score, scaled_score)
    } else {
      lst[[name]] <- value  # If it's not a data frame, just store it as-is
    }
  }
  lst
}


#' Get top regulons by absolute deltaPagoda value
#'
#' Selects the top \`n\` regulons for a given condition based on the
#' largest absolute \`DeltaPagoda\` values, then returns them ordered
#' by their signed \`DeltaPagoda\`.
#'
#' @param data A data frame containing at least the columns
#'   \`"Comparison"\`, \`"DeltaPagoda"\`, and \`"TF"\`.
#' @param condition Character. The condition name to filter on
#'   (matches the \`Comparison\` column).
#' @param top_n Integer. The number of top regulons (by absolute
#'   \`DeltaPagoda\`) to return.
#'
#' @return A character vector of regulon names (\`TF\`), length
#'   \`<= top_n\`, ordered by ascending \`DeltaPagoda\`.
#'
#'
#' @importFrom dplyr filter arrange slice_head pull
#' @export
get_top_tfs <- function(data, condition, top_n) {
  data %>%
    filter(Comparison == condition, !is.na(DeltaPagoda)) %>%
    arrange(desc(abs(DeltaPagoda))) %>%
    slice_head(n = top_n) %>%
    arrange(DeltaPagoda) %>%
    pull(TF)
}


#' Create a pseudo-logarithmic transformation
#'
#' Constructs a bidirectional transformation that behaves like a logarithm
#' but handles zero and negative values gracefully via a signed log1p.
#' Useful for plotting scales that span negative and positive values.
#'
#' @param base Numeric. The logarithm base. Defaults to 10.
#'
#' @return A \`trans\` object (from the **scales** package) with
#'   \`transform\`, \`inverse\`, and \`domain\` defined for signed data.
#'
#'
#' @importFrom scales trans_new
#' @export
pseudo_log_trans <- function(base = 10) {
  trans_new(
    name = paste0("pseudo_log", base),
    transform = function(x) sign(x) * log1p(abs(x)) / log(base),
    inverse = function(x) sign(x) * (base^abs(x) - 1),
    domain = c(-Inf, Inf)
  )
}


#' Calculate Percentage Change
#'
#' Vectorised helper that returns the percentage change from a baseline
#' (`old_val`) to an updated value (`new_val`).
#'
#' The function handles division-by-zero explicitly:
#' * If `old_val` is `0` and `new_val` is positive, the result is `Inf`.
#' * If both `old_val` and `new_val` are `0`, the result is `0`.
#'
#' @param new_val Numeric vector of new or current values.
#' @param old_val Numeric vector of baseline values (same length as `new_val`).
#'
#' @return Numeric vector of percentage changes. Values can be finite,
#'   `Inf`, or `0`, matching the length of the inputs.
#'
#' @examples
#' calculate_pct_change(c(120, 80, 0), c(100, 100, 0))
#' #> 20 -20 Inf
#'
#' @export
calculate_pct_change <- function(new_val, old_val) {
  change <- ifelse(old_val == 0,
                   ifelse(new_val > 0, Inf, 0),
                   ((new_val - old_val) / old_val) * 100)
  return(change)
}


#' Initialize Empty Data-Frame Containers
#'
#' Creates and returns a named list of two empty data.frames, \code{edges} and
#' \code{combined_data}, for downstream wrangling or visualization.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{edges}}{An empty \code{data.frame} to store edge lists.}
#'   \item{\code{combined_data}}{An empty \code{data.frame} to accumulate joined data.}
#' }
#' @export
#'
#' @examples
#' df_list <- initialize_data_frames()
#' str(df_list)
initialize_data_frames <- function() {
  list(
    edges = data.frame(),
    combined_data = data.frame()
  )
}

#' Generate Divergent Color Vector for log₂FC Values
#'
#' Given a numeric vector of log₂ fold-change values, this returns a character
#' vector of the same length with colors interpolated between
#' \code{"cornflowerblue"} → \code{"white"} → \code{"coral1"}. Missing or
#' sentinel values (\code{-999}) are mapped to \code{"white"}.
#'
#' @param log2fc_values Numeric vector of log₂ fold-change values (may contain NAs or -999).
#' @param num_colors    Integer; number of discrete colors to generate (default 100).
#'
#' @return Character vector of colors, one per element of \code{log2fc_values}.
#' @export
generate_log2fc_colors <- function(log2fc_values, num_colors = 100) {
  valid_log2fc <- log2fc_values[log2fc_values != -999 & !is.na(log2fc_values)]
  max_abs_log2fc <- max(abs(valid_log2fc), na.rm = TRUE)
  breaks <- seq(-max_abs_log2fc, max_abs_log2fc, length.out = num_colors + 1)
  color_palette <- colorRampPalette(c("cornflowerblue", "white", "coral1"))(num_colors)

  colors <- cut(log2fc_values, breaks = breaks, labels = color_palette, include.lowest = TRUE)
  colors <- as.character(colors)
  colors[log2fc_values == -999 | is.na(log2fc_values)] <- "white"

  colors
}


#' Set Vertex Aesthetics for an igraph Object
#'
#' Assigns fill‐colors, sizes, and label text sizes to the vertices of a graph
#' based on precomputed log₂FC color mappings and a set of “core” regulon nodes.
#'
#' @param g               An \code{igraph} graph object whose vertices are named.
#' @param log2fc_colors   Named character vector of colors (names = vertex names)
#'                        to apply based on log₂ fold-change values.
#' @param regulons        Character vector of vertex names considered “core” regulons.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{colors}}{Character vector of vertex fill‐colors.}
#'   \item{\code{sizes}}{Numeric vector of vertex sizes (e.g.\ radius).}
#'   \item{\code{label_cex}}{Numeric vector of label text scaling factors.}
#' }
#' @export
set_vertex_attributes <- function(g, log2fc_colors, regulons) {
  vertex_colors <- rep("white", vcount(g))
  vertex_colors[V(g)$name %in% names(log2fc_colors)] <- log2fc_colors[V(g)$name %in% names(log2fc_colors)]
  vertex_colors[V(g)$name %in% regulons] <- "darkgoldenrod1"

  vertex_size <- ifelse(V(g)$name %in% regulons, 10, 5)
  vertex_label_cex <- ifelse(V(g)$name %in% regulons, 1.0, 0.6)

  list(colors = vertex_colors, sizes = vertex_size, label_cex = vertex_label_cex)
}

#' Retrieve PubMed Article Count for a Gene
#'
#' Queries NCBI PubMed via \pkg{rentrez} and returns the total number of
#' articles matching the given gene symbol or search term.
#'
#' @param gene Character scalar: gene symbol or PubMed search term.
#' @return Integer: total count of PubMed records for that term.
#' @importFrom rentrez entrez_search
#' @export
#'
#' @examples
get_n_pubmed_articles_per_gene <- function(gene) {
    search_result <- entrez_search(db = "pubmed", term = gene)
    return(search_result$count)
  }