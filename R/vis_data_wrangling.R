#' Count the Number of Outliers in a Dataframe
#'
#' This function calculates the number of outliers in the 'prioritization_score'
#' column of a dataframe. An outlier is defined as a value that is more than a
#' certain number of standard deviations away from the mean. This function counts
#' how many values exceed one, two, and three standard deviations above the mean.
#'
#' @param df A dataframe that must contain a column named 'prioritization_score'.
#'
#' @return A vector with three elements, each representing the count of values
#' that are more than one, two, and three standard deviations above the mean, respectively.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'your_dataframe' has a column 'prioritization_score'
#'   outlier_counts <- getNumberOfOutliers(your_dataframe)
#'   print(outlier_counts)
#' }
#'
#' @export
#'
#' @importFrom stats mean sd
#' @importFrom base sum colnames stop
getNumberOfOutliers <- function(df) {
  # Checking if 'prioritization_score' column exists
  if(!"prioritization_score" %in% colnames(df)) {
    stop("The dataframe does not contain a 'prioritization_score' column.")
  }

  # Calculating required statistics
  mean_val <- mean(df$prioritization_score)
  sd_val <- sd(df$prioritization_score)
  threshold_1 <- mean_val + sd_val
  threshold_2 <- mean_val + 2 * sd_val
  threshold_3 <- mean_val + 3 * sd_val


  # Counting elements above each threshold
  count_above_threshold_1 <- sum(df$prioritization_score > threshold_1)
  count_above_threshold_2 <- sum(df$prioritization_score > threshold_2)
  count_above_threshold_3 <- sum(df$prioritization_score > threshold_3)

  return(c(count_above_threshold_1,count_above_threshold_2,count_above_threshold_3))
}

#' Scale Prioritization Score
#'
#' This function scales the prioritization scores in a dataframe column.
#' If all values are non-negative, they are scaled between 0 and 1. If negative
#' values are present, they are scaled to a range between -1 and 1.
#'
#' @param df A dataframe containing the prioritization scores to be scaled.
#' @param score_column A string specifying the name of the column that contains
#' the prioritization scores.
#'
#' @return A dataframe with an additional column `scaled_score` containing the
#' scaled prioritization scores.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'your_dataframe' has a score column named 'priority_score'
#'   scaled_df <- scale_prioritization_score(your_dataframe, 'priority_score')
#'   print(scaled_df)
#' }
#'
#' @export
#'
#' @importFrom stats min max
#' @import dplyr
#' @importFrom dplyr mutate
#' @importFrom base any is.na ifelse
scale_prioritization_score <- function(df, score_column) {
  # Check if any value in the specified score column is negative
  has_negatives <- any(df[[score_column]] < 0, na.rm = TRUE)

  # Determine the scaling method based on the presence of negative values
  scaled_score <- if (!has_negatives) {
    # Scale between 0 and 1 for non-negative values
    (df[[score_column]] - min(df[[score_column]])) /
      (max(df[[score_column]]) - min(df[[score_column]]))
  } else {
    # Scale between -1 and 1 for ranges that include negative values
    ifelse(
      df[[score_column]] >= 0,
      df[[score_column]] / max(df[[score_column]]),
      df[[score_column]] / abs(min(df[[score_column]]))
    )
  }

  # Add the scaled score to the data frame
  df <- df %>% mutate(scaled_score = scaled_score)

  return(df)
}

#' Fill Missing Ligand Data for Cytosig Comparison
#'
#' This function identifies ligands present in `matching_genes_lset` but missing
#' in `method_results`. It then creates new rows with these ligands, setting
#' their interaction score to zero, receptor to "unknown", and the sender cluster
#' to "mixed". These new rows are appended to the original `method_results`
#' dataframe.
#'
#' @param method_results A dataframe containing the results of a method with
#' columns for ligand, interaction, score, and receiver_cluster.
#' @param matching_genes_lset A vector of ligand gene symbols to compare against
#' the method results.
#'
#' @return A dataframe that combines the original `method_results` with new rows
#' for missing ligands, ensuring that all ligands in `matching_genes_lset` are
#' represented in the output.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'cytosig_results' is your method results dataframe
#'   # and 'ligand_set' is your set of ligand gene symbols
#'   completed_results <- fillGapsForCytosigComparison(cytosig_results, ligand_set)
#'   print(completed_results)
#' }
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#' @importFrom base setdiff unique data.frame paste0
fillGapsForCytosigComparison <- function(method_results,matching_genes_lset){

  # Vector of new ligands
  new_ligands <- setdiff(matching_genes_lset, method_results$ligand)

  # Generate new rows for the missing ligands
  new_rows <- data.frame(
    interaction = paste0(new_ligands, "-unknown"),
    score = 0,
    ligand = new_ligands,
    receptor = "unknown",
    sender_cluster = "mixed",
    receiver_cluster = unique(method_results$receiver_cluster)  # This assumes all rows have the same receiver cluster
  )

  # Combine the original data frame with the new rows
  result_df <- bind_rows(method_results, new_rows)
  return(result_df)
}

#' Count Overlapping Items Based on Rank Threshold
#'
#' Computes the number of items in a dataframe that have both 'rank.x' and 'rank.y'
#' less than a specified threshold. This function is useful for identifying the
#' count of items that meet a certain criteria of overlap between two ranking
#' systems within a dataset.
#'
#' @param i A numeric value specifying the threshold for the rank. Items with
#' 'rank.x' and 'rank.y' both less than this value are counted as overlapping.
#' @param data A dataframe that contains at least two columns: 'rank.x' and 'rank.y',
#' which represent ranking systems or scores for comparison.
#'
#' @return An integer representing the number of items that have both 'rank.x'
#' and 'rank.y' less than the specified threshold `i`.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'your_dataframe' contains columns 'rank.x' and 'rank.y'
#'   overlap_count <- count_overlaps(10, your_dataframe)
#'   print(overlap_count)
#' }
#'
#' @export
count_overlaps <- function(i, data) {
  sum(data$rank.x < i & data$rank.y < i)
}


#' Generate Comparison Object from NicheNet Analysis
#'
#' This function processes the results from a NicheNet analysis, stored within
#' a list of tibbles, to format them for comparison purposes. It renames columns
#' to standardize ligand and score nomenclature, and adds new columns for receiver
#' cluster, receptor (set to "unknown"), interaction (combining ligand and receptor),
#' and sender cluster (set to "mixed").
#'
#' @param nichenetObject A list containing NicheNet analysis results, where each
#' element is a tibble with at least the columns `test_ligand` and `aupr_corrected`.
#'
#' @return A modified version of the input list where each tibble has been processed
#' to include 'ligand', 'score', 'receiver_cluster', 'receptor', 'interaction', and
#' 'sender_cluster' columns, with 'receptor' and 'sender_cluster' set to default values.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'nichenet_results' is your list of NicheNet analysis results
#'   comparison_object <- generateComparisonObjectFromNicheNet(nichenet_results)
#'   print(comparison_object)
#' }
#'
#' @export
#' @importFrom dplyr rename select mutate
#' @importFrom stats setNames
generateComparisonObjectFromNicheNet <- function(nichenetObject){

  if (is.null(nichenetObject)) return(NULL)

  # Load the required library
  library(dplyr)

  # Modify column names for each tibble inside each list
  nichenet_results_modified <- lapply(names(nichenet_results), function(name) {
    df <- nichenet_results[[name]]

    df %>%
      rename(
        ligand = test_ligand,
        score = aupr_corrected
      ) %>%
      select(
        ligand,
        score
      ) %>%
      mutate(
        receiver_cluster = name,
        receptor = "unknown",
        interaction = paste(ligand, receptor, sep = "-"),
        sender_cluster = "mixed"
      )
  })

  # Update the original nichenet_results with the modified version
  return(nichenet_results_modified)

}

#' Rename Decipher Score Column
#'
#' This function renames the 'decipher_score' column in a given dataframe to 'score'.
#' It's particularly useful for standardizing column names across different dataframes
#' for consistent processing or analysis.
#'
#' @param df A dataframe that contains a column named 'decipher_score'.
#'
#' @return A dataframe identical to the input but with the 'decipher_score' column
#' renamed to 'score'.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'your_dataframe' contains a column named 'decipher_score'
#'   updated_df <- renameDecipherScore(your_dataframe)
#'   print(updated_df)
#' }
#'
#' @export
#' @importFrom dplyr rename
renameDecipherScore <- function(df) {
  df <- df %>%
    rename(score = decipher_score)
  return(df)
}


#' Preprocess NATMI Dataframe
#'
#' This function takes a dataframe containing NATMI analysis results and performs
#' a series of preprocessing steps including filtering, renaming, and reordering of
#' columns, as well as adding a new interaction column. It concludes by scaling
#' the prioritization scores for comparability.
#'
#' @param natmi_df A dataframe containing NATMI analysis results. The dataframe
#' must include columns for ligand symbols, receptor symbols, sending and target
#' clusters, and log2 transformed fold change of edge expression weight.
#'
#' @return A dataframe where infinite values in the log2 transformed fold change
#' of edge expression weight have been removed, column names have been standardized,
#' an interaction column has been added, and only relevant columns are retained.
#' Furthermore, the prioritization score is scaled for comparability, and the dataframe
#' is ordered by receiver and prioritization score in descending order.
#'
#' @examples
#' \dontrun{
#'   processed_natmi_df <- preProcessNATMI(natmi_results_all)
#'   print(processed_natmi_df)
#' }
#'
#' @export
#' @import dplyr
#' @importFrom stats desc
preProcessNATMI <- function(natmi_df){
  #NULL check
  if (is.null(natmi_df)) return(NULL)

  result <- natmi_df %>%
    #filter
    filter(Log2.transformed.fold.change.of.edge.expression.weight != Inf,
           Log2.transformed.fold.change.of.edge.expression.weight!= -Inf) %>%
    #rename
    rename(ligand = Ligand.symbol,
           receptor = Receptor.symbol,
           sender = Sending.cluster,
           receiver = Target.cluster,
           prioritization_score = Log2.transformed.fold.change.of.edge.expression.weight) %>%
    #add interaction column
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    #select columns of interest for downstream analysis
    select(sender, receiver, interaction, prioritization_score) %>%
    #organized based on score
    arrange(receiver, desc(prioritization_score))

  #scale prioritization score to make it comparable
  result <- scale_prioritization_score(result,"prioritization_score")
  return(result)
}

#' Preprocess Connectome Dataframe
#'
#' This function preprocesses a dataframe containing connectome data. It filters out
#' entries with infinite scores, renames columns for standardization, adds an
#' interaction column concatenating ligand and receptor names, selects relevant columns,
#' and arranges the dataframe based on the receiver and descending prioritization scores.
#' Finally, it scales the prioritization scores to make them comparable.
#'
#' @param connectome_df A dataframe representing connectome data, expected to include
#' 'source', 'target', and 'score' columns, with 'source' and 'target' representing
#' the sender and receiver, respectively, and 'score' indicating the interaction strength.
#'
#' @return A processed dataframe with infinite scores filtered out, columns renamed to
#' 'sender', 'receiver', and 'prioritization_score', an additional 'interaction' column
#' created, and data arranged by receiver and prioritization score in descending order.
#' The 'prioritization_score' is scaled for comparability.
#'
#' @examples
#' \dontrun{
#'   processed_connectome_df <- preProcessConnectome(connectome_data)
#'   print(processed_connectome_df)
#' }
#'
#' @export
#' @import dplyr
#' @importFrom stats desc
preProcessConnectome <- function(connectome_df){
  #NULL check
  if (is.null(connectome_df)) return(NULL)

  result <- connectome_df%>%
    #filter
    filter(score!=Inf) %>%
    #rename
    rename(sender = source,
           receiver = target,
           prioritization_score = score) %>%
    #add interaction column
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    #select columns of interest for downstream analysis
    select(sender, receiver, interaction, prioritization_score) %>%
    #organized based on score
    arrange(receiver, desc(prioritization_score))

  #scale prioritization score to make it comparable
  result <- scale_prioritization_score(result,"prioritization_score")

  return(result)
}

#' Preprocess LIANA Dataframe
#'
#' Processes a dataframe resulting from LIANA (Ligand-Receptor Inference And Network
#' Analysis) by adding an interaction column that concatenates ligand and receptor
#' names, renaming columns for clarity, selecting specific columns for downstream
#' analysis, and arranging the entries based on the receiver and the prioritization
#' score in descending order. The prioritization scores are then scaled for
#' comparability.
#'
#' @param liana_df A dataframe containing ligand-receptor interaction data from LIANA,
#' including columns for 'ligand', 'receptor', 'source' (sender), 'target' (receiver),
#' and 'interaction_stat' (the statistical score representing the strength or significance
#' of the interaction).
#'
#' @return A processed dataframe with a new 'interaction' column, standardized column
#' names ('sender', 'receiver', 'prioritization_score'), and entries arranged by receiver
#' and prioritization score in descending order. The 'prioritization_score' is scaled
#' for comparability across different datasets or analysis methods.
#'
#' @examples
#' \dontrun{
#'   processed_liana_df <- preProcessLIANA(liana_data)
#'   print(processed_liana_df)
#' }
#'
#' @export
#' @importFrom dplyr mutate rename select arrange
#' @importFrom stats desc
preProcessLIANA <- function(liana_df){

  #NULL check
  if (is.null(liana_df)) return(NULL)

  result <-  liana_df %>%
    # add interaction column
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    #rename columns
    dplyr::rename(
      sender = source,
      receiver = target,
      prioritization_score = interaction_stat
    )  %>%
    #select columns of interest
    select(sender, receiver, interaction, prioritization_score) %>%
    #order by score
    arrange(receiver, desc(prioritization_score))

  #scale score
  result <- scale_prioritization_score(result,"prioritization_score")
  return(result)

}


preProcessDecipher <- function(decipher_results) {
  if (is.null(decipher_results)) return(NULL)

  decipher_bound <- bind_rows(decipher_results)

  result <- decipher_bound %>%
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    rename(
      receiver = receiver_cluster,
      sender = sender_cluster,
      prioritization_score = decipher_score
    )

  result <- scale_prioritization_score(result, "prioritization_score")
  return(result)
}

preProcessNicheNet <- function(nichenet_prior_table_all_clusters) {
  if (is.null(nichenet_prior_table_all_clusters)) return(NULL)

  nichenet_bound <- bind_rows(nichenet_prior_table_all_clusters)

  result <- nichenet_bound %>%
    mutate(interaction = paste(ligand, receptor, sep = "-"))

  result <- scale_prioritization_score(result, "prioritization_score")
  return(result)
}




#' Prepare Data for Correlation Analysis
#'
#' Groups data by interaction and receiver, then filters to retain only the rows where
#' the absolute value of the prioritization score is equal to the maximum absolute
#' prioritization score within each group. This process identifies the most significant
#' interactions for each receiver-interaction pair. The dataset is then ungrouped,
#' and only relevant columns are retained. Finally, the data is ordered by prioritization
#' score in ascending order.
#'
#' @param df A dataframe that must contain at least three columns: `interaction`,
#' `receiver`, and `prioritization_score`. The `interaction` column represents the
#' unique interactions between ligands and receptors. The `receiver` column specifies
#' the entities receiving the interaction, and `prioritization_score` quantifies the
#' significance or strength of each interaction.
#'
#' @return A dataframe filtered to include only the most significant interaction for
#' each interaction-receiver pair, with rows ordered by prioritization score. The
#' resulting dataframe contains three columns: `interaction`, `receiver`, and
#' `prioritization_score`.
#'
#' @examples
#' \dontrun{
#'   significant_interactions <- prepareDataForCorrelationAnalysis(your_dataframe)
#'   print(significant_interactions)
#' }
#'
#' @export
#' @importFrom dplyr group_by filter ungroup select arrange
prepareDataForCorrelationAnalysis <- function(df){
  result <- df %>%
    group_by(interaction, receiver) %>%
    filter(abs(prioritization_score) == max(abs(prioritization_score))) %>%
    ungroup()%>%
    select(interaction,receiver,prioritization_score) %>%
    arrange(prioritization_score)
  return(result)
}

#' Calculate Spearman Correlation and Search Space Between Interaction Methods
#'
#' This function computes the Spearman correlation coefficient for the prioritization
#' scores between all unique pairs of interaction methods provided. It also determines
#' the search space (the minimum number of top-ranked interactions) needed to achieve
#' at least 100 overlaps between the pairs of methods. The results are returned in two
#' matrices: one for the Spearman correlation coefficients and another for the search
#' space sizes.
#'
#' @param method_results_list A named list of dataframes, each representing the
#' results from a different interaction method. Each dataframe should have columns
#' named `interaction`, `receiver`, and `prioritization_score`, with `prioritization_score.x`
#' and `prioritization_score.y` being used for merged dataframes to represent scores
#' from the two methods being compared.
#'
#' @return A list containing two matrices: `spearman` with Spearman correlation coefficients
#' for all pairs of methods, and `k_matrix` with the corresponding search space size
#' needed to achieve at least 100 overlapping interactions between each pair of methods.
#'
#' @examples
#' \dontrun{
#'   method_list <- list(method1 = df1, method2 = df2)
#'   correlation_and_search_space <- getInteractionCorrelationAndSearchSpaceBetweenMethods(method_list)
#'   print(correlation_and_search_space$spearman)
#'   print(correlation_and_search_space$k_matrix)
#' }
#'
#' @export
#' @importFrom dplyr inner_join mutate
#' @importFrom stats cor rank
#' @importFrom base which list matrix names print paste sapply
getInteractionCorrelationAndSearchSpaceBetweenMethods <- function(method_results_list){

  spearman_matrix <- matrix(NA, nrow = length(method_results_list),
                            ncol = length(method_results_list),
                            dimnames = list(names(method_results_list), names(method_results_list)))

  k_matrix <- matrix(NA, nrow = length(method_results_list),
                     ncol = length(method_results_list),
                     dimnames = list(names(method_results_list), names(method_results_list)))

  # Iterate over all unique pairs of data frames
  for (name1 in names(method_results_list)) {
    for (name2 in names(method_results_list)) {
      #if (name1 < name2) {  # This ensures each pair is only considered once
      print(paste("Comparing", name1, "and", name2))

      # Inner join the data frames
      merged_data <- inner_join(method_results_list[[name1]], method_results_list[[name2]],
                                by = c("interaction", "receiver")) %>%
        mutate(rank.x = rank(-prioritization_score.x, ties.method = "min"),
               rank.y = rank(-prioritization_score.y, ties.method = "min"))

      overlap_counts <- sapply(100:max(c(merged_data$rank.x, merged_data$rank.y)),
                               "count_overlaps",
                               data = merged_data)

      # Find the smallest i that satisfies the condition
      smallest_i_index <- which(overlap_counts >= 100)[1]
      smallest_i <- smallest_i_index + 99  # Adjusting the index to match the actual i value

      # Subset for calculating Spearman correlation
      ind <- which(merged_data$rank.x < smallest_i & merged_data$rank.y < smallest_i)

      # Calculate Spearman correlation
      spearman_cor <- cor(merged_data$rank.x[ind], merged_data$rank.y[ind], method = "spearman")

      print(paste("Spearman correlation for", name1, "and", name2, ":", spearman_cor))
      # Update the matrix
      spearman_matrix[name1, name2] <- spearman_cor
      k_matrix[name1,name2] <- smallest_i
      #}
    }
  }
  result <- list()
  result[["spearman"]] <- spearman_matrix
  result[["k_matrix"]] <- k_matrix
  return(result)
}


getoverlapTable <- function(method_results_list){
  result <- data.frame(
    index = 1,
    comparison = dataset_name,
    reported_interactions_decipher = nrow(method_results_list$Decipher),
    reported_interactions_liana = nrow(method_results_list$`LIANA+`),
    reported_interactions_nichenet = nrow(method_results_list$NicheNet),
    reported_interactions_connectome = nrow(method_results_list$Connectome),
    reported_interactions_natmi = nrow(method_results_list$NATMI),
    above_1_sd_decipher = getNumberOfOutliers(method_results_list$Decipher)[1],
    above_1_sd_liana = getNumberOfOutliers(method_results_list$`LIANA+`)[1],
    above_1_sd_nichenet = getNumberOfOutliers(method_results_list$NicheNet)[1],
    above_1_sd_connectome = getNumberOfOutliers(method_results_list$Connectome)[1],
    above_1_sd_natmi = getNumberOfOutliers(method_results_list$NATMI)[1],
    above_2_sd_decipher = getNumberOfOutliers(method_results_list$Decipher)[2],
    above_2_sd_liana = getNumberOfOutliers(method_results_list$`LIANA+`)[2],
    above_2_sd_nichenet = getNumberOfOutliers(method_results_list$NicheNet)[2],
    above_2_sd_connectome = getNumberOfOutliers(method_results_list$Connectome)[2],
    above_2_sd_natmi = getNumberOfOutliers(method_results_list$NATMI)[2],
    above_3_sd_decipher = getNumberOfOutliers(method_results_list$Decipher)[3],
    above_3_sd_liana = getNumberOfOutliers(method_results_list$`LIANA+`)[3],
    above_3_sd_nichenet = getNumberOfOutliers(method_results_list$NicheNet)[3],
    above_3_sd_connectome = getNumberOfOutliers(method_results_list$Connectome)[3],
    above_3_sd_natmi = getNumberOfOutliers(method_results_list$NATMI)[3]
  )
  return(result)
}

#' Extract Top N Unique Interactions
#'
#' This function sorts a given dataframe of cell-cell interactions (CCIs) by the
#' absolute value of their prioritization scores in descending order. It then extracts
#' the unique interaction identifiers (from the `interaction` column) for the top N
#' interactions based on this sorting.
#'
#' @param cci_df A dataframe containing cell-cell interactions. It must include a
#' column named `prioritization_score` for scoring the interactions and a column named
#' `interaction` containing interaction identifiers.
#' @param n_top An integer specifying the number of top interactions to extract after
#' sorting by prioritization score.
#'
#' @return A vector of interaction identifiers for the top N unique interactions, based
#' on the absolute values of their prioritization scores.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'cci_dataframe' is your dataframe containing CCIs with columns
#'   # 'prioritization_score' and 'interaction'
#'   top_interactions <- getSet(cci_dataframe, 100)
#'   print(top_interactions)
#' }
#'
#' @export
#' @importFrom dplyr arrange pull unique
getSet <- function(cci_df,n_top){
  result <- cci_df %>%
    arrange(-abs(prioritization_score)) %>%
    pull(interaction) %>%
    unique()
  result <- result[1:n_top]
  return(result)
}


getInteractionOverlap <- function(method_results_list){
  interaction_counts <- table(unlist(method_results_list))

  # Convert to a data frame
  interaction_counts_df <- as.data.frame(interaction_counts)
  # Rename the columns
  names(interaction_counts_df) <- c("Interaction", "Count")
  # View the data frame

  meta_overlap <- list()
  for(this_method in names(method_results_list)){
    this_method_interactions <- method_results_list[[this_method]]
    ind_df <- match(this_method_interactions, interaction_counts_df$Interaction)
    counts <- interaction_counts_df$Count[ind_df]
    method_overlap <- table(counts)
    method_overlap_df <- as.data.frame(method_overlap)
    # Rename the columns
    names(method_overlap_df) <- c("overlap", "Count")
    method_overlap_df$method <- this_method
    meta_overlap[[this_method]] <- method_overlap_df
  }
  meta_overlap <- bind_rows(meta_overlap)
  return(meta_overlap)
}


#' Summarize Z-Scores Across Clusters
#'
#' This function calculates the median Z-score for each ligand across different
#' clusters or samples, provided a list of Z-score files. Each file is expected
#' to contain Z-scores for all ligands in a particular cluster or sample.
#'
#' @param z_score_files A vector of filenames (including paths) pointing to the CSV
#' files containing Z-scores for ligands. Each file represents a different cluster
#' or sample.
#' @param z_score_folder folder where the CSV files containing Z-scores are located
#' @param mapping_table a mapping data frame to map from cytosig ligand proteins to gene names
#'
#' @return A list where each element is named after the input file (minus the '.csv'
#' extension) and contains a vector of median Z-scores for each ligand within that
#' file.
#'
#' @examples
#' \dontrun{
#'   z_score_files <- c("/path/to/cluster1_z_scores.csv", "/path/to/cluster2_z_scores.csv")
#'   median_z_scores <- summarizeZScores(z_score_files)
#'   print(median_z_scores)
#' }
#'
#' @importFrom data.table fread
#' @importFrom stats median apply
#' @importFrom utils file.path
#' @importFrom stringr str_remove
#' @importFrom tibble enframe
#' @importFrom tidyr pivot_wider drop_na
#' @importFrom dplyr mutate left_join
summarizeZScores <- function(z_score_files,z_score_folder,mapping_table){
  z.median <- list()
  for(this_cluster in c(1:length(z_score_files))){
    z_score_filename <- z_score_files[this_cluster]
    filename <- stringr::str_remove(z_score_filename,pattern = ".csv")

    sample_z <- fread(file.path(z_score_folder,z_score_filename),header = TRUE)
    ligand_names <- sample_z$V1
    sample_z <- sample_z[,-1]
    sample_z <- as.matrix(sample_z)
    rownames(sample_z) <- ligand_names

    z.median[[filename]] <- apply(sample_z, 1, median)
  }

  cytosig_significance <- z.median %>%
    tibble::enframe(name = "rowname", value = "values") %>%
    tidyr::unnest_longer(values, indices_to = "ligand") %>%
    tidyr::pivot_wider(names_from = rowname, values_from = values) %>%
    left_join(mapping_table,by = c("ligand"="ligand")) %>%
    mutate(gene = if_else(gene == "",NA,gene))%>%
    tidyr::drop_na(gene)

  return(cytosig_significance)
}


#' Prepare LIANA Data for Cytosig Comparison
#'
#' Transforms a LIANA data frame for comparison with Cytosig results. This involves
#' selecting relevant columns, renaming them for consistency with Cytosig data
#' standards, and creating an interaction column by concatenating ligand and receptor
#' names. The data is then split into a list of data frames based on the receiver
#' cluster, facilitating further comparison or analysis.
#'
#' @param liana_df A dataframe from LIANA analysis containing columns for ligand,
#' receptor, source cluster, target cluster, and interaction statistic.
#'
#' @return A list of data frames, each corresponding to a different receiver cluster
#' from the original `liana_df`. Columns are standardized for direct comparison with
#' Cytosig results, and an additional 'interaction' column is included, combining
#' ligand and receptor names.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'liana_df' is your dataframe obtained from LIANA analysis
#'   liana_prepared <- prepareLianaForCytosigComparison(liana_df)
#'   # Now 'liana_prepared' is ready for comparison with Cytosig results
#' }
#'
#' @importFrom dplyr select rename mutate split
#' @importFrom base paste return
prepareLianaForCytosigComparison <- function(liana_df){
  if (is.null(liana_df)) return(NULL)

  liana_df <-liana_df %>%
    select(ligand,receptor,source,target,interaction_stat)%>%
    rename(
      sender_cluster = source,
      receiver_cluster = target,
      score = interaction_stat
    ) %>%
    mutate(interaction = paste(ligand,receptor,sep = "-"))
  liana_results_for_comparison <- split(liana_df, liana_df$receiver_cluster)
  return(liana_results_for_comparison)
}

#' Prepare Connectome Data for Cytosig Comparison
#'
#' Processes a dataframe from Connectome analysis by filtering out infinite scores,
#' selecting essential columns, renaming columns for consistency with Cytosig format,
#' and creating an 'interaction' column by concatenating ligand and receptor names.
#' The dataframe is then split into a list of dataframes based on the receiver cluster,
#' making it compatible for comparison or analysis alongside Cytosig results.
#'
#' @param connectome_df A dataframe containing interaction data from Connectome analysis,
#' including ligand, receptor, source cluster, target cluster, and interaction score.
#'
#' @return A list of dataframes, each corresponding to a different receiver cluster
#' in the `connectome_df`. Columns are aligned with Cytosig standards, and an 'interaction'
#' column is added that combines ligand and receptor names for easy comparison.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'connectome_df' is your dataframe from Connectome analysis
#'   connectome_prepared <- prepareConnectomeForCytosigComparison(connectome_df)
#'   # 'connectome_prepared' is now structured for comparison with Cytosig results
#' }
#'
#' @importFrom dplyr filter select rename mutate split
#' @importFrom base paste return
prepareConnectomeForCytosigComparison <- function(connectome_df){
  if (is.null(connectome_df)) return(NULL)
  connectome_df <-connectome_df %>%
    filter(score!=Inf & score!=-Inf)%>%
    select(ligand,receptor,source,target,score)%>%
    rename(
      sender_cluster = source,
      receiver_cluster = target
    ) %>%
    mutate(interaction = paste(ligand,receptor,sep = "-"))
  connectome_results_for_comparison <- split(connectome_df, connectome_df$receiver_cluster)
  return(connectome_results_for_comparison)
}

#' Prepare NATMI Data for Cytosig Comparison
#'
#' Filters and reformats a NATMI analysis dataframe for compatibility with Cytosig comparison.
#' This involves filtering out infinite values from the log2 transformed fold change of edge
#' expression weights, renaming columns to standard ligand-receptor interaction terms, creating
#' an interaction identifier by concatenating ligand and receptor names, and finally splitting
#' the dataframe into a list of dataframes based on the receiver cluster.
#'
#' @param natmi_df A dataframe resulting from NATMI analysis containing columns for ligand symbols,
#' receptor symbols, source and target clusters, and log2 transformed fold change of edge
#' expression weights.
#'
#' @return A list of dataframes, each keyed by the receiver cluster, formatted for direct comparison
#' with Cytosig results. This includes standardized column names (`ligand`, `receptor`, `sender_cluster`,
#' `receiver_cluster`, `score`) and a new `interaction` column.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'natmi_df' is your dataframe from NATMI analysis
#'   natmi_prepared <- prepareNatmiForCytosigComparison(natmi_df)
#'   # 'natmi_prepared' is now ready for comparison with Cytosig results
#' }
#'
#' @importFrom dplyr filter rename mutate split
#' @importFrom base paste return
prepareNatmiForCytosigComparison <- function(natmi_df){
  if (is.null(natmi_df)) return(NULL)
  natmi_df <- natmi_df %>%
    filter(Log2.transformed.fold.change.of.edge.expression.weight != Inf,
           Log2.transformed.fold.change.of.edge.expression.weight!= -Inf) %>%
    rename(ligand = Ligand.symbol,
           receptor = Receptor.symbol) %>%
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    rename(sender_cluster = Sending.cluster,
           receiver_cluster = Target.cluster,
           score = Log2.transformed.fold.change.of.edge.expression.weight)
  natmi_results_for_comparison <- split(natmi_df, natmi_df$receiver_cluster)
  return(natmi_results_for_comparison)
}


#' Extract Predictions and Responses for Comparison Across Methods
#'
#' This function compares predictions from various methods against response data
#' available in a significance dataset. It processes each method's predictions, aligns
#' them with response data based on matching genes, and then generates plots to
#' visualize the comparison. The function aims to assess how closely the predictions
#' of each method match actual response data, considering only significant genes and
#' those within a specified ligand set.
#'
#' @param results_to_compare A named list of dataframes, each containing the prediction
#' results of a different method. Each dataframe should include columns for ligand
#' names and scores.
#' @param cytosig_significance A dataframe containing significance data for genes,
#' used as the basis for response data in the comparison.
#' @param L.set A dataframe specifying the set of ligands considered in the analysis,
#' used to filter the predictions and responses.
#' @param seurat_object_oi A Seurat object or equivalent containing gene expression
#' data, used to further filter the genes considered in the analysis based on their
#' presence in the dataset.
#'
#' @return A list containing two elements: `predictions` and `responses`, each a list
#' where each element corresponds to the prediction vectors or response vectors for
#' a method, respectively.
#'
#' @examples
#' \dontrun{
#'   # Example usage
#'   result <- getPredictionsResponsesForMethods(results_to_compare,
#'                                               cytosig_significance,
#'                                               L.set,
#'                                               seurat_object_oi)
#' }
#'
#' @importFrom dplyr filter group_by ungroup pull
#' @importFrom stats setNames
#' @importFrom graphics png dev.off plot
#' @importFrom utils file.path paste
#' @importFrom base c list match names return
getPredictionsResponsesForMethods <- function(results_to_compare, cytosig_significance,L.set,seurat_object_oi,output_figures_filepath){
  all_predictions_across_methods <- list()
  all_responses_across_methods <- list()

  #TODO: fails when one of the methods is missing a cell-type, figure out a way to add a filler!

  for(this_method in names(results_to_compare)){
    comparison_method <- this_method
    comparison_method_results <- results_to_compare[[this_method]]

    #ROC
    all_predictions <- c()
    all_responses <- c()
    for(this_receiver_ct in names(results_to_compare$Decipher)){
      if(this_receiver_ct %in% names(comparison_method_results)){
        method_results <- comparison_method_results[[this_receiver_ct]]

        matching_ligands <- intersect(method_results$ligand,cytosig_significance$gene)
        prediction_vector <- method_results$score
        names(prediction_vector) <- method_results$ligand
        #first I need to remove genes that don't exist in our seurat object - these are not fair
        matching_genes <- intersect(rownames(seurat_object_oi),cytosig_significance$gene)
        #then I need to remove genes that are not in our L set, these are also not fair
        matching_genes_lset <- intersect(L.set$ligand,matching_genes)
        method_results_for_cytosig <- fillGapsForCytosigComparison(method_results,matching_genes_lset)
        duplicated(method_results_for_cytosig$ligand)

        #keep comparison ligands
        method_results_for_cytosig <- method_results_for_cytosig %>% filter(ligand %in% matching_genes_lset)

        # keep single ligands
        method_results_for_cytosig_filtered <- method_results_for_cytosig %>%
          group_by(ligand) %>%
          filter(abs(score) == max(abs(score))) %>%
          ungroup()

        prediction_vector <- method_results_for_cytosig_filtered %>% pull(score)
        response_vector <- cytosig_significance %>% filter(gene %in% matching_genes_lset) %>% pull(!!this_receiver_ct)

        # Extract the ligand labels from method_results_for_cytosig_filtered
        prediction_vector_labels <- method_results_for_cytosig_filtered %>% pull(ligand)

        # Extract the gene labels from cytosig.results
        response_vector_labels <- cytosig_significance %>%
          filter(gene %in% matching_genes_lset) %>%
          pull(gene)

        # Determine the common ordering of labels
        prediction_vector_new_order <- match(response_vector_labels, prediction_vector_labels)

        # Reorder both vectors based on the common ordering
        prediction_vector_ordered <- prediction_vector[prediction_vector_new_order]

        all_predictions <- c(all_predictions, prediction_vector_ordered)
        all_responses <- c(all_responses, response_vector)

        max_x <- max(abs(response_vector))
        max_y <- max(abs(prediction_vector_ordered))
        # png(file.path(output_figures_filepath,paste(comparison_method,this_receiver_ct,"comparison.test.png")))
        # p <- plot(response_vector,prediction_vector_ordered,xlim = c(-1*max_x,max_x),ylim = c(-1*max_y,max_y),main=this_receiver_ct)
        # print(p)
        # dev.off()


      } else {
        method_results <- comparison_method_results[[1]]
        method_results$score <- 0
        method_results$receiver_cluster <- this_receiver_ct

        matching_ligands <- intersect(method_results$ligand,cytosig_significance$gene)
        prediction_vector <- method_results$score
        names(prediction_vector) <- method_results$ligand
        #first I need to remove genes that don't exist in our seurat object - these are not fair
        matching_genes <- intersect(rownames(seurat_object_oi),cytosig_significance$gene)
        #then I need to remove genes that are not in our L set, these are also not fair
        matching_genes_lset <- intersect(L.set$ligand,matching_genes)
        method_results_for_cytosig <- fillGapsForCytosigComparison(method_results,matching_genes_lset)
        duplicated(method_results_for_cytosig$ligand)

        #keep comparison ligands
        method_results_for_cytosig <- method_results_for_cytosig %>% filter(ligand %in% matching_genes_lset)

        # keep single ligands
        method_results_for_cytosig_filtered <- method_results_for_cytosig %>%
          group_by(ligand) %>%
          filter(abs(score) == max(abs(score))) %>%
          ungroup()

        prediction_vector <- method_results_for_cytosig_filtered %>% pull(score)
        response_vector <- cytosig_significance %>% filter(gene %in% matching_genes_lset) %>% pull(!!this_receiver_ct)

        # Extract the ligand labels from method_results_for_cytosig_filtered
        prediction_vector_labels <- method_results_for_cytosig_filtered %>% pull(ligand)

        # Extract the gene labels from cytosig.results
        response_vector_labels <- cytosig_significance %>%
          filter(gene %in% matching_genes_lset) %>%
          pull(gene)

        # Determine the common ordering of labels
        prediction_vector_new_order <- match(response_vector_labels, prediction_vector_labels)

        # Reorder both vectors based on the common ordering
        prediction_vector_ordered <- prediction_vector[prediction_vector_new_order]

        all_predictions <- c(all_predictions, prediction_vector_ordered)
        all_responses <- c(all_responses, response_vector)

        max_x <- max(abs(response_vector))
        max_y <- max(abs(prediction_vector_ordered))
        png(file.path(output_figures_filepath,paste(comparison_method,this_receiver_ct,"comparison.test.png")))
        p <- plot(response_vector,prediction_vector_ordered,xlim = c(-1*max_x,max_x),ylim = c(-1*max_y,max_y),main=this_receiver_ct)
        print(p)
        dev.off()
      }

    }

    all_predictions_across_methods[[this_method]] <- all_predictions
    all_responses_across_methods[[this_method]] <- all_responses

  }
  result <- list()
  result[["predictions"]] <- all_predictions_across_methods
  result[["responses"]] <- all_responses_across_methods
  return(result)
}

#' Plot ROC Curves and Extract AUC for Multiple Methods
#'
#' For each specified threshold, this function plots the Receiver Operating
#' Characteristic (ROC) curve and calculates the Area Under the Curve (AUC) for
#' predictions versus responses across multiple methods. The function iteratively
#' generates and saves a ROC curve plot for each threshold and accumulates the
#' AUC values for each method and threshold combination.
#'
#' @param predictions A list of prediction vectors, named by method, where each
#' element represents predicted scores for a set of instances.
#' @param responses A list of response vectors, named by method, where each element
#' contains the actual outcomes (e.g., binary response data) corresponding to
#' the predictions. The structure should match that of `predictions`.
#' @param output_figures_filepath A character string specifying the path where ROC
#' curve plots should be saved.
#'
#' @return A list of AUC values for each method and threshold combination. Each
#' element of the list is named by the threshold value and contains a nested list
#' with AUC values for each method.
#'
#' @examples
#' \dontrun{
#'   predictions = list(Method1 = c(0.1, 0.4, 0.35, 0.8),
#'                      Method2 = c(0.2, 0.3, 0.5, 0.9))
#'   responses = list(Method1 = c(0, 0, 1, 1),
#'                    Method2 = c(0, 1, 0, 1))
#'   output_path = "path/to/save/plots/"
#'   auc_values = plotROCAndExtractAUC(predictions, responses, output_path)
#'   print(auc_values)
#' }
#'
#' @importFrom graphics png dev.off plot lines text
#' @importFrom utils file.path paste
#' @importFrom pROC roc auc
plotROCAndExtractAUC <- function(predictions,responses,output_figures_filepath,dataset_name){
  # Calculate the max threshold
  #max_threshold <- floor(max(abs(responses$Decipher)))
  max_threshold <- floor(max(responses$Decipher))
  # Generate a sequence of values
  seq_values <- seq(0, max_threshold, by = 1)
  auc_for_meta <- list()

  for(this_threshold in seq_values){
    png(file.path(output_figures_filepath,paste(this_threshold,"roc_comparison.png",sep="_")),width = 15,height=15,units="cm",res = 400)
    # Initialize the plot
    plot(1,
         type="n",
         xlab="FPR",
         ylab="TPR",
         xlim=c(1, 0),
         ylim=c(0, 1),
         main=paste(dataset_name," (",">",this_threshold,"\u03c3)",sep=""),
         cex.lab = 1.5)

    # Initialize colors
    colors <- c("coral1", "chartreuse4", "cornflowerblue","hotpink2","mediumorchid3")

    # Plot ROC curves
    i <- 1
    auc_methods <- list()
    for (method in names(predictions)) {
      # Calculate binary responses based on threshold
      binary_responses <- as.numeric(responses[[method]] > this_threshold)

      # Calculate ROC curve
      roc_curve <- roc(binary_responses, predictions[[method]],quiet=TRUE)

      # Plot ROC curve
      lines(roc_curve, col=colors[i], lwd=2)

      # Add text for AUC
      auc_value <- round(auc(roc_curve), 3)
      n_true <- sum(binary_responses)
      #text(0.85, 1-(i-1)*0.1, paste(method, ": AUC = ", auc_value, sep=""), col=colors[i])
      text(0.3, 0.4-(i-1)*0.1, paste(method, " AUC = ", auc_value, sep=""), col=colors[i],cex=1.1)

      i <- i + 1

      random_auc <- data.frame(
        threshold = character(),
        auc = numeric()
      )

      #auc_methods[[method]] <- auc_value
      auc_methods[[method]] <- list(auc = auc_value, n_true = n_true)
    }
    auc_for_meta[[as.character(this_threshold)]] <- auc_methods
    dev.off()
  }
  return(auc_for_meta)
}


#' Get Top Regulons Based on Delta Threshold
#'
#' Filters and identifies the top regulons from a dataframe based on a specified
#' delta threshold. Regulons with their `deltaPagoda` absolute values exceeding this
#' threshold are selected and ordered by `deltaPagoda` in descending order. This function
#' is useful for identifying the most significant regulons by their changes in Pagoda
#' scores, possibly indicating a significant regulatory effect or change.
#'
#' @param df A dataframe containing regulon data, expected to include columns
#' `deltaPagoda` for regulon scores and `name` for regulon names.
#' @param delta_threshold A numeric value specifying the minimum absolute value of
#' `deltaPagoda` required for a regulon to be considered as "top" or significant.
#'
#' @return A vector of unique regulon names that meet the specified delta threshold,
#' ordered by their `deltaPagoda` scores in descending order.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'regulon_df' is your dataframe and you're interested in regulons
#'   # with a deltaPagoda score above 0.5
#'   top_regulons <- get_top_regulons(regulon_df, 0.5)
#'   print(top_regulons)
#' }
#'
#' @importFrom dplyr filter arrange pull unique
#' @export
get_top_regulons <- function(df,delta_threshold) {
  top_regulons <- df %>%
    filter(abs(deltaPagoda) > delta_threshold) %>%
    arrange(-deltaPagoda) %>%
    pull(name) %>%
    unique()
  return(top_regulons)
}

#' Extract Top Regulons and Corresponding Cell Types
#'
#' This function processes a list of dataframes containing regulon delta values by cluster,
#' filters for "real" entries, identifies the top regulons based on a specified delta threshold,
#' and extracts the corresponding cell types. It returns a list containing the unique top regulons
#' across all clusters and the names of clusters (cell types) associated with these top regulons.
#'
#' @param regulon_deltas_by_cluster A list of dataframes, each representing regulon
#' delta values for a specific cell type or cluster. Each dataframe is expected to have
#' been structured to be compatible with `filter_real` and `get_top_regulons` functions.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{top_regulons}: A vector of unique top regulon names across all clusters,
#'   filtered based on a delta threshold.
#'   \item \code{cts}: The names of clusters (cell types) corresponding to the identified
#'   top regulons.
#' }
#'
#' @examples
#' \dontrun{
#'   regulon_deltas <- list(
#'     Cluster1 = data.frame(deltaPagoda = c(4, 2, 1), name = c("RegulonA", "RegulonB", "RegulonC")),
#'     Cluster2 = data.frame(deltaPagoda = c(5, 1, 3), name = c("RegulonA", "RegulonD", "RegulonE"))
#'   )
#'   result <- getTopRegulonsAndCts(regulon_deltas)
#'   print(result$top_regulons)
#'   print(result$cts)
#' }
#'
#' @importFrom stats unlist unique
#' @seealso \code{\link{filter_real}}, \code{\link{get_top_regulons}}
getTopRegulonsAndCts <- function(regulon_deltas_by_cluster){
  filtered_data <- lapply(regulon_deltas_by_cluster, filter_real)
  filtered_data_top <- lapply(filtered_data, get_top_regulons,delta_threshold=3)
  top_regulons <- unique(unlist(filtered_data_top))
  # Get the list of cell types
  cts_corresponding_to_top_regulons <- names(filtered_data_top)
  result <- list()
  result[["top_regulons"]] <- top_regulons
  result[["cts"]] <- cts_corresponding_to_top_regulons
  return(result)
}

#' Create a Matrix of Regulon Delta Values Across Cell Types
#'
#' Constructs a matrix where each row represents a top regulon and each column
#' represents a cell type (CT). The values in the matrix are the delta Pagoda
#' scores for the corresponding regulon in each CT. If a regulon is not present
#' in a CT, its delta Pagoda score in the matrix is set to 0.
#'
#' @param regulon_deltas_all_clusters A list where each element is a dataframe
#' containing the delta Pagoda scores for regulons within a specific cell type.
#' Each list element's name corresponds to a cell type, and each dataframe must
#' contain at least the columns `name` for the regulon names and `deltaPagoda`
#' for their scores.
#' @param top_regulons_and_cts A list with two elements: `top_regulons`, a vector
#' containing the names of the top regulons to include in the matrix, and `cts`,
#' a vector containing the names of the cell types to be included as columns in
#' the matrix.
#'
#' @return A matrix with rows corresponding to top regulons and columns to cell types.
#' The entries of the matrix are the delta Pagoda scores of the regulons in each cell
#' type, with 0s where a regulon is not present in a cell type.
#'
#' @examples
#' \dontrun{
#'   regulon_deltas <- list(
#'     CT1 = data.frame(name = c("Regulon1", "Regulon2"), deltaPagoda = c(2.3, -1.5)),
#'     CT2 = data.frame(name = c("Regulon2", "Regulon3"), deltaPagoda = c(3.1, 0.5))
#'   )
#'   top_regulons_and_cts <- list(
#'     top_regulons = c("Regulon1", "Regulon2", "Regulon3"),
#'     cts = c("CT1", "CT2")
#'   )
#'   delta_matrix <- pull_top_regulons_cts_delta_matrix(regulon_deltas, top_regulons_and_cts)
#'   print(delta_matrix)
#' }
#'
#' @export
pull_top_regulons_cts_delta_matrix <- function(regulon_deltas_all_clusters,top_regulons_and_cts){
  result_matrix <- matrix(0,
                          nrow = length(top_regulons_and_cts$top_regulons),
                          ncol = length(top_regulons_and_cts$cts),
                          dimnames = list(top_regulons_and_cts$top_regulons, top_regulons_and_cts$cts))

  for(this_ct in colnames(result_matrix)){
    regulon_deltas_this_ct <- regulon_deltas_all_clusters[[this_ct]]
    for(regulon in top_regulons_and_cts$top_regulons){
      if(regulon %in% regulon_deltas_this_ct$name){
        result_matrix[regulon,this_ct] <- regulon_deltas_this_ct$deltaPagoda[which(rownames(regulon_deltas_this_ct) == regulon)]
      }
    }
  }

  return(result_matrix)
}



#' Get Differential Expression Markers for Ligands and Receptors by Cell Type
#'
#' This function aggregates differential expression markers for ligands and receptors
#' across different cell types specified in the input list. It constructs a data frame
#' that includes average log fold change, the cluster or cell type, the gene name,
#' and the percentage of cells in two conditions expressing the gene.
#'
#' @param lr.marker.list A list where each element represents a different cell type
#' and contains a data frame of ligand-receptor differential expression markers for that cell type.
#' Each data frame should have row names set to gene names and contain columns for
#' average log2 fold change (avg_log2FC), and percentage of cells expressing the gene in
#' two different conditions (pct.1, pct.2).
#'
#' @return A data frame consolidating all input differential expression markers with columns:
#' \itemize{
#' \item \code{avg_log2FC}: Average log2 fold change of expression.
#' \item \code{cluster}: Cell type or cluster.
#' \item \code{gene}: Gene name.
#' \item \code{pct.1}: Percentage of cells in condition 1 expressing the gene.
#' \item \code{pct.2}: Percentage of cells in condition 2 expressing the gene.
#' }
#'
#' @examples
#' # Assuming 'lr_marker_list' is pre-defined with the appropriate structure
#' results <- getLigandReceptorDiffExprMarkersByCt(lr_marker_list)
#'
#' @importFrom yourPackageName addEntryToDF
#' @export
getLigandReceptorDiffExprMarkersByCt <- function(lr.marker.list){
  condition_lr_markers <- data.frame(
    avg_log2FC = 0,
    cluster = "",
    gene = "",
    pct.1 = 0,
    pct.2 = 0
  )
  for(this_cluster in names(lr.marker.list)){
    this_lr_markers <- lr.marker.list[[this_cluster]]
    this_lr_markers$cluster <- this_cluster
    this_lr_markers$gene <- rownames(this_lr_markers)
    #ind.ligand <- which(this_lr_markers$gene == "ADAM17")
    this_cond_lr_markers <- this_lr_markers[,c("avg_log2FC","cluster","gene","pct.1","pct.2")]
    condition_lr_markers <- addEntryToDF(condition_lr_markers,this_cond_lr_markers)
  }
  condition_lr_markers <- condition_lr_markers[-1,]
  return(condition_lr_markers)
}

#' Calculate P-Value for a Given Value Against a Set of Base Values
#'
#' This function calculates the two-sided p-value of a given value based on its
#' standard score (z-score) in relation to a set of base values. It assumes that the
#' base values are normally distributed. The function computes the mean and standard
#' deviation of the base values, calculates the z-score of the given value, and then
#' derives the p-value using the standard normal distribution.
#'
#' @param base_values A numeric vector of base values assumed to be from a normal distribution.
#' @param real_value A single numeric value for which the p-value is to be calculated.
#'
#' @return A numeric value representing the two-sided p-value.
#'
#' @examples
#' base_values <- rnorm(100, mean = 50, sd = 10)
#' real_value <- 60
#' p_value <- calculate_p_value(base_values, real_value)
#' print(p_value)
#'
#' @export
calculate_p_value <- function(base_values,real_value){
  mu <- mean(base_values)
  sigma <- sd(base_values)
  z <- (real_value - mu)/sigma
  p_value <- p_value <- 2 * (1 - pnorm(abs(z)))
  return(p_value)
}


#' Perform t-tests by feature, grouped by a factor
#'
#' This function conducts t-tests for each feature in a data matrix, comparing values across groups defined by a grouping factor.
#'
#' @param data_matrix A matrix where rows represent samples and columns represent features.
#' @param group_factor A factor specifying the grouping of samples.
#'
#' @return A named vector of p-values, where each p-value corresponds to a feature in the data matrix.
#'
#' @details This function iterates over each feature in the data matrix, conducts a t-test comparing the values across groups defined by the grouping factor, and returns the p-values for each feature.
#'
#' @examples
#' data_matrix <- matrix(rnorm(100), nrow = 10)
#' group_factor <- factor(rep(c("A", "B"), each = 5))
#' do_t_test_by_feature_by_grouping_factor(data_matrix, group_factor)
#'
#' @importFrom stats t.test
#'
#' @export
do_t_test_by_feature_by_grouping_factor <- function(data_matrix,group_factor){
  # Initialize a vector to store the p-values
  p_values <- numeric(nrow(data_matrix))

  # Loop over each feature
  for (i in 1:nrow(data_matrix)) {
    # Extract the values for this feature
    feature_values <- data_matrix[i, ]

    # Perform a t-test comparing the case and control values
    t_test_result <- t.test(feature_values ~ group_factor)

    # Store the p-value for this test
    p_values[i] <- t_test_result$p.value
  }

  # If you want to preserve the feature names in the result
  feature_names <- rownames(data_matrix)
  p_values_named <- setNames(p_values, feature_names)
  return(p_values_named)
}

#' Get Top Overlapping Regulons
#'
#' This function processes a list of regulon delta values by cluster, applies filtering, and identifies the
#' top overlapping regulons based on their frequency across clusters. It filters real data points, identifies top regulons
#' based on a delta threshold, and then counts the frequency of each regulon across clusters to determine
#' the most common ones.
#'
#' @param regulon_deltas_by_cluster A list where each element corresponds to a cluster and contains a data frame of
#'        regulons and their delta values. These data frames are expected to have a structure that can be processed
#'        by the `filter_real` and `get_top_regulons` functions.
#' @param n_overlapping The number of top overlapping regulons to return.
#'
#' @return A vector of the top `n_overlapping` regulon names, ordered by the frequency of their appearance across
#'         clusters, from the most to the least frequent.
#'
#' @examples
#' # Assume regulonDeltas is a list of data frames with delta values for different clusters
#' topRegulons <- getTopOverlappingRegulons(regulonDeltas, 5)
#' print(topRegulons)
#'
#' @export
getTopOverlappingRegulons <- function(regulon_deltas_by_cluster,n_overlapping){
  filtered_data <- lapply(regulon_deltas_by_cluster, filter_real)
  filtered_data_top <- lapply(filtered_data, get_top_regulons,delta_threshold=3)
  all_tfs <- unlist(filtered_data_top)
  tf_counts <- table(all_tfs)
  top_regulons <- names(tf_counts[order(tf_counts,decreasing=TRUE)])[1:n_overlapping]
  return(top_regulons)
}

#' Calculate Spearman Correlation Matrix for Ligand-Receptor-Transcription Factor Axes
#'
#' This function computes a Spearman correlation matrix for decipher scores by regulon and cluster,
#' optionally filtered by top regulons. It creates a matrix showing the correlation of signed weights
#' calculated as the product of the sign of Spearman correlation coefficients and a permuted importance measure
#' across different clusters for given regulon combinations defined by ligand-receptor pairs.
#'
#' @param decipher_scores_by_regulon_and_cluster A list of data frames where each list element corresponds
#'        to a cluster, and each data frame contains columns for ligands, receptors, regulons,
#'        Spearman correlation coefficients ('spearman.cor'), and permuted importance ('imp.perm').
#' @param top_regulons An optional vector of regulon names to focus the analysis on specific regulons.
#'        If NULL, the analysis will include all regulons in the data.
#'
#' @return A list containing two matrices:
#'         - `$spearman_matrix`: A matrix of Spearman correlation coefficients.
#'         - `$label_matrix`: A matrix of counts of unique regulons compared between clusters.
#'
#' @examples
#' # Example usage:
#' decipher_scores <- list(
#'   Cluster1 = data.frame(ligand = c("L1", "L2"), receptor = c("R1", "R2"), regulon = c("Reg1", "Reg2"),
#'                         spearman.cor = c(0.5, -0.5), imp.perm = c(1.2, 0.8)),
#'   Cluster2 = data.frame(ligand = c("L1", "L3"), receptor = c("R1", "R3"), regulon = c("Reg1", "Reg3"),
#'                         spearman.cor = c(0.6, -0.6), imp.perm = c(1.0, 0.9))
#' )
#' result <- GetSpearmanLRTF(decipher_scores)
#' print(result$spearman_matrix)
#' print(result$label_matrix)
#'
#' @importFrom dplyr mutate filter inner_join
#' @importFrom tidyr paste
#' @export
GetSpearmanLRTF <- function(decipher_scores_by_regulon_and_cluster,top_regulons = NULL){
  blank_matrix <- matrix(NA, nrow = length(decipher_scores_by_regulon_and_cluster),
                         ncol = length(decipher_scores_by_regulon_and_cluster),
                         dimnames = list(names(decipher_scores_by_regulon_and_cluster), names(decipher_scores_by_regulon_and_cluster)))
  spearman_matrix <- blank_matrix
  label_matrix <- blank_matrix

  for(this_ct in names(decipher_scores_by_regulon_and_cluster)){
    for(this_ct_2 in names(decipher_scores_by_regulon_and_cluster)){
      CT1 <- decipher_scores_by_regulon_and_cluster[[this_ct]]
      CT2 <- decipher_scores_by_regulon_and_cluster[[this_ct_2]]

      CT1 <- CT1  %>% mutate(lr_tf_axis = paste(ligand, receptor, regulon, sep = '-'),
                             signed_weight = sign(spearman.cor)*imp.perm)
      CT2 <- CT2  %>% mutate(lr_tf_axis = paste(ligand, receptor, regulon, sep = '-'),
                             signed_weight = sign(spearman.cor)*imp.perm)

      if(!is.null(top_regulons)){
        CT1 <- CT1 %>% filter(regulon %in% top_regulons) %>% mutate(lr_tf_axis = paste(ligand, receptor, regulon, sep = '-'),
                                                                    signed_weight = sign(spearman.cor)*imp.perm)
        CT2 <- CT2 %>% filter(regulon %in% top_regulons) %>% mutate(lr_tf_axis = paste(ligand, receptor, regulon, sep = '-'),
                                                                    signed_weight = sign(spearman.cor)*imp.perm)
      }

      to_compare <- inner_join(CT1,CT2,by="lr_tf_axis")
      label_matrix[this_ct,this_ct_2] <- length(unique(to_compare$regulon.y))
      this_cor <- cor(to_compare$signed_weight.x,to_compare$signed_weight.y,method = "spearman")
      spearman_matrix[this_ct,this_ct_2] <- this_cor
    }
  }

  spearman_matrix[is.na(spearman_matrix)] <- 0

  result <- list()
  result[["spearman_matrix"]] <- spearman_matrix
  result[["label_matrix"]] <- label_matrix
  return(result)
}


#' Select Top Rows Based on Method and Score
#'
#' This function selects the top rows from a dataset based on the specified method and scoring criteria.
#' It can operate in either a positive or negative mode, determined by the `negative` flag.
#' In positive mode, it selects the top positive scores. In negative mode, and if a specific method is
#' selected, it picks the top negative scores for that method.
#'
#' @param method_results Data frame containing method results, including columns for 'receiver', 'method',
#'        and 'prioritization_score'.
#' @param n_rows The number of top rows to select based on the score (default is 1).
#' @param negative Logical flag indicating whether to select for negative scores (default is FALSE, select
#'        for positive scores).
#' @param selected_method Optionally specify a method to filter by in negative selection mode.
#'
#' @return A data frame of the selected top rows with an additional column 'label' that indicates
#'         whether the selection was for the top positive or top negative scores.
#'
#' @examples
#' # Example DataFrame
#' df <- data.frame(
#'   receiver = c("A", "B", "A", "B"),
#'   method = c("Method1", "Method1", "Method2", "Method2"),
#'   prioritization_score = c(0.9, 0.85, 0.95, 0.80)
#' )
#' # Get top positive rows by default
#' top_positive <- getTopRowsByMethod(df)
#' # Get top negative rows for a specific method
#' top_negative <- getTopRowsByMethod(df, negative = TRUE, selected_method = "Method2")
#'
#' @importFrom dplyr filter group_by slice_max slice_min mutate ungroup
#' @export
getTopRowsByMethod <- function(method_results,n_rows=1,negative = FALSE,selected_method=NULL){
  if(!negative){
    selected_rows_positive <- bound_rows %>%
      group_by(receiver,method)  %>%
      slice_max(n = 1, order_by = prioritization_score) %>%
      mutate(label = paste("Top positive",method)) %>%
      ungroup()
  } else if (negative){
    if(!is.null(selected_method)){
      selected_rows_negative <- bound_rows %>%
        filter(method == selected_method) %>%
        group_by(receiver,method)  %>%
        slice_min(n = 1, order_by = prioritization_score)%>%
        mutate(label = paste("Top negative",method)) %>%
        ungroup()
    }

  }
}

#' Construct Active Signaling Networks Using NicheNet
#'
#' This function uses NicheNet to identify active signaling networks within a Seurat object under
#' specific conditions. It filters ligands, targets, and transcription factors to construct a network
#' based on differential expression analysis between specified conditions.
#'
#' @param seurat_object_oi A Seurat object containing the data for analysis.
#' @param parameters A list containing parameters such as the condition of interest (`condition_oi`)
#'        and the reference condition (`condition_reference`).
#' @param df_for_comparison A data frame with columns 'receiver' and 'ligand' specifying which
#'        ligand-receiver pairs to analyze.
#' @param lr_network A data frame representing the ligand-receptor network with columns 'from' and 'to'.
#' @param ligand_target_matrix A matrix where rows are genes and columns are potential ligands, used
#'        to identify target genes of ligands.
#' @param ligand_tf_matrix A matrix where rows are genes and columns are transcription factors, used
#'        to identify target genes of transcription factors.
#' @param weighted_networks A list containing weighted networks, specifically a 'lr_sig' data frame and
#'        potentially others like 'gr' for gene regulation.
#'
#' @return A list of active signaling networks, each corresponding to the ligand-receiver pairs analyzed.
#'
#' @examples
#' # Assuming you have a prepared Seurat object and necessary matrices and data frames:
#' results <- getNicheNetNetwork(seurat_object_oi, params, df_pairs, lr_network,
#'                               ligand_target_matrix, ligand_tf_matrix, networks)
#'
#' @importFrom dplyr distinct inner_join filter rownames_to_column
#' @importFrom Seurat Idents DefaultAssay SetIdent FindMarkers
#' @importFrom tibble rownames_to_column
#' @importFrom nichenetr get_ligand_signaling_path
#' @importFrom utils setTxtProgressBar txtProgressBar close
#' @export
getNicheNetNetwork <- function(seurat_object_oi, parameters,df_for_comparison,lr_network,ligand_target_matrix,ligand_tf_matrix,weighted_networks){
  Idents(seurat_object_oi) <- seurat_object_oi$cluster
  table(seurat_object_oi$condition)
  condition_oi = parameters$condition_oi
  condition_reference = parameters$condition_reference
  DefaultAssay(seurat_object_oi) <- "RNA"

  # For newer Seurat versions, you may need to run the following
  seurat_object_oi <- UpdateSeuratObject(seurat_object_oi)
  lr_network = lr_network %>% distinct(from, to)
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))


  # Define the total number of iterations
  total <- dim(df_for_comparison)[1]

  # Set up the progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  active_signaling_network_list <- list()
  for(i in 1:total){

    # Update the progress bar
    setTxtProgressBar(pb, i)

    this_receiver <- df_for_comparison$receiver[i]
    this_ligand <- df_for_comparison$ligand[i]

    seurat_obj_receiver= subset(seurat_object_oi, idents = this_receiver)
    seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["condition", drop=TRUE]])

    condition_oi = condition_oi
    condition_reference = condition_reference

    DE_table_receiver = Seurat::FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% tibble::rownames_to_column("gene")

    geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    geneset_oi = geneset_oi[which(geneset_oi %in% rownames(ligand_tf_matrix))]
    geneset_oi <- geneset_oi %>% .[. %in% weighted_networks$gr$to]

    ligands_all = this_ligand
    targets_all = geneset_oi

    active_signaling_network = nichenetr::get_ligand_signaling_path(
      ligand_tf_matrix = ligand_tf_matrix,
      ligands_all = ligands_all,
      targets_all = targets_all,
      weighted_networks = weighted_networks)

    active_signaling_network_list[[i]] <- active_signaling_network
  }

  # Close the progress bar
  return(active_signaling_network_list)
  close(pb)
}

#' Merge Selected Interactions from NicheNet Analysis
#'
#' This function merges results from NicheNet signaling network analysis with extended comparison data.
#' It integrates signaling data based on ligand-receiver interactions, annotating them with additional
#' information from an extended comparison data frame. The function processes each row individually,
#' extracting specific interaction data, and merges additional attributes like prioritization scores and labels.
#'
#' @param df_for_comparison A data frame specifying ligand-receiver pairs to analyze.
#' @param active_signaling_network_list A list containing NicheNet signaling network results, indexed by
#'        concatenation of ligand and receiver names.
#' @param df_for_comparison_extended An extended version of the `df_for_comparison` data frame that includes
#'        additional columns such as prioritization scores, method, rank, and label.
#'
#' @return A data frame with merged results, including columns for interaction, regulon, weight,
#'         receiver, method, regulon valuation, and possibly other attributes inherited from the extended comparison data frame.
#'
#' @examples
#' # Assuming df_for_comparison, active_signaling_network_list, and df_for_comparison_extended are pre-defined:
#' merged_results <- MergeSelectedInteractionsNicheNetNetworkNicheNetResults(
#'   df_for_comparison,
#'   active_signaling_network_list,
#'   df_for_comparison_extended
#' )
#'
#' @importFrom dplyr rowwise mutate filter transmute ungroup left_join
#' @importFrom tidyr unnest
#' @export
MergeSelectedInteractionsNicheNetNetworkNicheNetResults <- function(df_for_comparison,active_signaling_network_list,df_for_comparison_extended){
  merged_nichenet_results <- df_for_comparison %>%
    rowwise() %>%
    mutate(
      interaction_data = list(active_signaling_network_list[[paste(ligand,receiver)]]$sig)
    )%>%
    tidyr::unnest(interaction_data) %>%
    dplyr::filter(ligand == from) %>%
    dplyr::transmute(
      interaction = interaction,
      regulon = to,
      weight = weight,
      receiver = receiver
    ) %>%
    mutate(method = "NicheNet",regulon.val = 0) %>%
    ungroup() %>%
    left_join(df_for_comparison_extended %>% select(receiver,interaction,prioritization_score,method,rank),by = c("receiver","interaction","method")) %>%
    left_join(df_for_comparison_extended %>% select(receiver,interaction,label),by = c("receiver","interaction"))
  return(merged_nichenet_results)

}

#' Prepare Decipher Results for Overlap Plot
#'
#' This function prepares and merges Decipher results with comparison data for plotting overlaps.
#' It iterates over each row of a comparison data frame, selects relevant data, joins it with
#' decipher scores filtered by the specific interaction, and then combines all the results.
#' The final output is structured for overlap plotting, specifically configured for visual comparison
#' of Decipher scores across different clusters.
#'
#' @param df_for_comparison A data frame containing the initial comparison data including
#'        method, receiver, interaction, ligand, receptor, prioritization score, label, and rank.
#' @param decipher_scores_by_regulon_and_cluster A list where each element corresponds to a cluster
#'        and contains a data frame of decipher scores for different regulons within that cluster.
#'
#' @return A data frame ready for plotting, including the columns method, interaction, regulon,
#'         weight, receiver, regulon valuation, prioritization score, label, and rank, all
#'         annotated with method 'Decipher'.
#'
#' @examples
#' # Assuming df_for_comparison and decipher_scores_by_regulon_and_cluster are predefined:
#' results <- prepareDecipherResultsForOverlapPlot(df_for_comparison, decipher_scores_by_regulon_and_cluster)
#' head(results)
#'
#' @importFrom dplyr select filter left_join bind_rows mutate
#' @export
prepareDecipherResultsForOverlapPlot <- function(df_for_comparison,decipher_scores_by_regulon_and_cluster){
  decipher_scores_by_cluster_for_2.6_comparisons_list <- list()
  for(this_row in 1:dim(df_for_comparison)[1]){
    this_data <- df_for_comparison[this_row,] %>% select(method,receiver,interaction,ligand,receptor,prioritization_score,label,rank)
    to_join <- decipher_scores_by_regulon_and_cluster[[this_data$receiver]] %>%
      filter(interaction == this_data$interaction) %>%
      select(interaction,regulon,imp.perm,spearman.cor,regulon.val)%>%
      mutate(weight=imp.perm*sign(spearman.cor))%>%
      select(interaction,regulon,weight,regulon.val)
    joined_data <- left_join(this_data,to_join,by = c("interaction"))
    decipher_scores_by_cluster_for_2.6_comparisons_list[[this_row]] <- joined_data
  }

  decipher_scores_by_cluster_for_2.6_comparisons <- bind_rows(decipher_scores_by_cluster_for_2.6_comparisons_list)
  decipher_scores_by_cluster_for_2.6_comparisons <- decipher_scores_by_cluster_for_2.6_comparisons  %>%
    mutate(method = "Decipher") %>%
    select(method,interaction,regulon,weight,receiver,regulon.val,prioritization_score,label,rank)
  return(decipher_scores_by_cluster_for_2.6_comparisons)

}



#' Prepare a data frame for correlation analysis
#'
#' If the input is `NULL`, returns `NULL`.  For the special case of the
#' “Decipher” method, selects and orders the pre-aggregated scores;
#' otherwise, delegates to `prepareDataForCorrelationAnalysis()`.
#'
#' @param name Character. The name of the method or dataset (e.g. `"Decipher"`).
#' @param df A data frame containing at least the columns required for
#'   correlation analysis.  May be `NULL`.
#'
#' @return A data frame ready for correlation analysis, or `NULL` if `df` was `NULL`.
#'
#'
#' @importFrom dplyr select arrange
#' @export
prepareForCorrelation <- function(name, df) {
  if (is.null(df)) return(NULL)
  
  if (name == "Decipher") {
    #different for Decipher since scores are already aggregated across sender cell types
    df %>%
      select(interaction, receiver, prioritization_score) %>%
      arrange(prioritization_score)
  } else {
    prepareDataForCorrelationAnalysis(df)
  }
}

#' Load and validate regulon deltaPagoda data
#'
#' Reads an RDS file containing a named list of data frames for various cell types,
#' filters it to only include the specified `cell_types`, and checks that each
#' element is a data frame with the required columns `"name"` and `"deltaPagoda"`.
#' Any errors during loading or validation emit a warning and return an empty list.
#'
#' @param file_path Character. Path to the `.rds` file containing the regulon data.
#' @param cell_types Character vector. Names of cell types to retain from the loaded list.
#'
#' @return A named list of data frames (one per cell type in `cell_types`), each
#'   with columns `"name"` and `"deltaPagoda"`. Returns an empty list if the file
#'   cannot be read or fails validation.
#'
#'
#' @export
load_regulon_data <- function(file_path, cell_types) {
  tryCatch({
      data_list <- readRDS(file_path)
      # Ensure it's filtered for selected cell types if the file contains more
      data_list <- data_list[intersect(names(data_list), cell_types)]
      # Add basic validation
      if(!is.list(data_list)) stop("Loaded data is not a list.")
      if(length(data_list) > 0) {
          first_el <- data_list[[1]]
          if(!is.data.frame(first_el) || !all(c("name", "deltaPagoda") %in% colnames(first_el))) {
              stop("Dataframe structure is incorrect. Needs 'name' and 'deltaPagoda' columns.")
          }
      }
      return(data_list)
  }, error = function(e) {
      warning(paste("Error loading or validating file:", file_path, "-", e$message))
      # Return an empty list or handle appropriately
      return(list())
  })
}

#' Retrieve a regulon’s deltaPagoda value for a given identity
#'
#' Safely extracts the \`deltaPagoda\` value for a specified \`regulon\` from
#' a named list of data frames (one per identity). Returns \`NA\` if the
#' identity or regulon is not present or if the data frame structure is invalid.
#'
#' @param identity A character or factor representing the identity (e.g., cell type).
#' @param regulon_list A named list where each element is a data frame containing
#'   columns \`"name"\` and \`"deltaPagoda"\`, indexed by identity names.
#' @param regulon Character. The regulon name whose \`deltaPagoda\` value should be retrieved.
#'
#' @return A numeric \`deltaPagoda\` value if found; otherwise \`NA\`.
#'
#'
#' @export
get_deltaPagoda <- function(identity, regulon_list, regulon) {
  identity_char <- as.character(identity)
  if (is.null(regulon_list) || !identity_char %in% names(regulon_list)) {
    # warning(paste("Identity", identity_char, "not found in provided regulon list. Returning NA."))
    return(NA)
  }
  df <- regulon_list[[identity_char]]
  if (is.null(df) || !is.data.frame(df) || !all(c("name", "deltaPagoda") %in% colnames(df))) {
    # warning(paste("Required columns ('name', 'deltaPagoda') missing or data invalid for identity", identity_char, ". Returning NA."))
    return(NA)
  }
  if (regulon %in% df$name) {
    # Handle potential multiple matches (shouldn't happen with unique names)
    return(df$deltaPagoda[df$name == regulon][1])
  } else {
    # warning(paste("Regulon", regulon, "not found for identity", identity_char, ". Returning NA."))
    return(NA)
  }
}

#' Find the absolute maximum deltaPagoda across conditions
#'
#' Iterates through a list of lists of data frames, each representing
#' deltaPagoda values for various conditions and cell types, and returns
#' the largest absolute deltaPagoda. If no valid values are found, returns 1.
#'
#' @param deltas_list_of_lists A list whose elements are themselves lists of
#'   data frames. Each data frame should contain a numeric column named
#'   \`deltaPagoda\`.
#'
#' @return A single numeric value: the maximum absolute
#'   \`deltaPagoda\` across all provided data frames, or \`1\` if none are valid.
#'
#'
#' @export
find_absolute_max <- function(deltas_list_of_lists) {
  max_val <- -Inf # Initialize with negative infinity
  for (cond_list in deltas_list_of_lists) {
      if (!is.null(cond_list) && length(cond_list) > 0) {
          cond_max <- sapply(cond_list, function(df) {
              if (!is.null(df) && is.data.frame(df) && "deltaPagoda" %in% colnames(df) && nrow(df) > 0) {
                  current_max <- max(abs(df$deltaPagoda), na.rm = TRUE)
                  # Handle case where all deltaPagoda are NA or df is empty after NA removal
                  if (is.infinite(current_max)) {
                    return(-Inf) # Return -Inf if no valid values
                  } else {
                    return(current_max)
                  }
              } else {
                  return(-Inf) # Return -Inf for invalid/empty dataframes
              }
          })
          # Filter out -Inf before taking the max for the condition
          valid_cond_max <- cond_max[is.finite(cond_max)]
          if (length(valid_cond_max) > 0) {
             max_val <- max(max_val, max(valid_cond_max, na.rm = TRUE), na.rm = TRUE)
          }
      }
  }
   # If max_val is still -Inf (no valid data found), return a default like 1 or NA
  return(ifelse(is.finite(max_val), max_val, 1))
}



#' Convert regulon deltas into long format
#'
#' Takes a named list of data frames (one per cluster) and a condition label,
#' filters each data frame to retain only rows where `class == "real"`, adds
#' `Cluster` and `Condition` columns, and binds them into a single data frame.
#'
#' @param regulon_deltas A named list of data frames.  Each element represents
#'   a cluster’s regulon data and must include a column `class`.
#' @param condition Character. A label (e.g. `"moderate"` or `"severe"`) to
#'   assign in the `Condition` column for all rows.
#'
#' @return A data frame in long format with:
#'   - All rows from each input data frame where `class == "real"`.  
#'   - A `Cluster` column giving the source list element name.  
#'   - A `Condition` column set to the provided `condition` value.
#'
#'
#' @importFrom dplyr bind_rows filter
#' @export
get_long_deltas <- function(regulon_deltas, condition) {
  bind_rows(lapply(names(regulon_deltas), function(ct) {
    df <- regulon_deltas[[ct]] %>%
      filter(class == "real")  # Keep only real TFs
    df$Cluster <- ct
    df$Condition <- condition
    return(df)
  }))
}




#' Extract the top N interaction identifiers from a method result
#'
#' Safely retrieves up to \`n\` interaction names from \`method_result\`,
#' which may be a data frame (or matrix/DataFrame) with columns
#' \`"interaction"\` and \`"prioritization_score"\`, or a character vector.
#' Returns an empty character vector if inputs are NULL, empty, or invalid.
#'
#' @param method_result A data frame, matrix, \`DataFrame\`, or character
#'   vector.  If data frame–like, must contain columns:
#'   \`"interaction"\` (identifiers) and \`"prioritization_score"\`
#'   (numeric scores).  If a character vector, its first \`n\` elements
#'   will be returned.
#' @param n Integer. The maximum number of interaction identifiers to return.
#'
#' @return A character vector of up to \`n\` interaction names sorted by
#'   descending \`prioritization_score\`.  Returns \`character(0)\` if
#'   \`method_result\` is NULL, empty, or improperly structured.
#'
#'
#' @export
getSet <- function(method_result, n) {
  # Basic Input Checks
  if (is.null(method_result)) {
    # warning("Input to getSet is NULL. Returning empty set.")
    return(character(0))
  }
  # Check if it's dataframe-like and has rows
  if (!is.data.frame(method_result) && !is(method_result, "DataFrame") && !is.matrix(method_result)) {
     # If it's some other structure, maybe it's already a vector of interactions?
     if(is.character(method_result)) {
         return(head(method_result, n))
     }
     warning("Input to getSet is not a recognized dataframe-like structure or character vector. Returning empty set.")
     return(character(0))
  }
   if (nrow(method_result) == 0) {
    # warning("Input to getSet has 0 rows. Returning empty set.")
    return(character(0))
  }

  # Check for 'interaction' column 
  interaction_col <- "interaction" 
  if (!interaction_col %in% colnames(method_result)) {
     warning(paste("Column '", interaction_col, "' not found in data for a method. Returning empty set.", sep=""))
     return(character(0))
  }

  # Extract top N unique interactions
  # Ensure the column is character type
  score_col <- "prioritization_score"
  sorted_indices <- order(method_result[[score_col]], decreasing = TRUE, na.last = TRUE)
  sorted_method_result <- method_result[sorted_indices, , drop = FALSE] # Use drop=FALSE for safety
  interactions <- unique(as.character(sorted_method_result[[interaction_col]]))
  return(head(interactions, n))
}


#' Calculate specific intersection counts across multiple method sets
#'
#' Given a named list of character vectors (interaction sets) for five methods
#' ("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome"), computes, for
#' every non-empty subset of these methods, the number of elements that are
#' present in all methods of the subset and absent from all methods not in the subset.
#'
#' @param list_input Named list of character vectors.  The names should be
#'   method identifiers (any or all of "Decipher", "NicheNet", "LIANA+",
#'   "NATMI", "Connectome"), and each element is a vector of interaction IDs.
#'
#' @return A data frame with two columns:
#'   - `Intersection_Name`: a string of method names separated by `" & "`, 
#'     representing the subset.  
#'   - `Count`: integer, the number of interactions found exclusively in that subset.
#'
#' @examples
#' \dontrun{
#' sets <- list(
#'   Decipher    = c("A", "B", "C"),
#'   NicheNet    = c("B", "C", "D"),
#'   `LIANA+`    = c("C", "E"),
#'   NATMI       = c("A", "C"),
#'   Connectome  = c("C", "F")
#' )
#' df_intersections <- calculate_all_intersections(sets)
#' # Look for the count of elements only in "Decipher & NATMI"
#' subset(df_intersections, Intersection_Name == "Decipher & NATMI")
#' }
#'
#' @export
calculate_all_intersections <- function(list_input) {
  method_names <- names(list_input)
  expected_methods <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")
  n_methods <- length(expected_methods)

  # Ensure all 5 methods are potentially present, use empty set if NULL or missing
  sets <- list()
  for(m in expected_methods) {
    sets[[m]] <- if (!is.null(list_input[[m]]) && length(list_input[[m]]) > 0) list_input[[m]] else character(0)
  }

  # List to store results: Intersection Name -> Count
  all_intersections_info <- list()

  # Loop through all possible intersection sizes (k = 1 to 5)
  for (k in 1:n_methods) {
    # Generate all combinations (subsets) of method names of size k
    combinations_k <- combn(expected_methods, k, simplify = FALSE)

    # Process each combination (subset S)
    for (subset_S in combinations_k) {
      # Sort subset for consistent naming
      subset_S_sorted <- sort(subset_S)
      # Create the intersection name
      intersection_name <- paste(subset_S_sorted, collapse = " & ")

      # Identify methods NOT in the current subset (subset NotS)
      subset_NotS <- setdiff(expected_methods, subset_S)

      # --- Calculate the count for elements ONLY in subset_S ---

      # 1. Find the intersection of all sets IN subset_S
      # Need to handle case where subset_S has only 1 element
      if (length(subset_S) == 1) {
        intersect_S <- sets[[subset_S[[1]]]]
      } else {
        intersect_S <- Reduce(intersect, sets[subset_S])
      }

      # If the intersection is already empty, the specific count is 0
      if (length(intersect_S) == 0) {
        count <- 0
      } else {
        # 2. Find the union of all sets NOT in subset_S (if any)
        union_NotS <- character(0) # Initialize as empty set
        if (length(subset_NotS) > 0) {
          union_NotS <- Reduce(union, sets[subset_NotS])
        }

        # 3. Find elements in intersect_S that are NOT in union_NotS
        specific_intersect_S_elements <- setdiff(intersect_S, union_NotS)
        count <- length(specific_intersect_S_elements)
      }

      # Store result
      all_intersections_info[[intersection_name]] <- count
    }
  }

  # Convert list to dataframe
  results_df <- data.frame(
    Intersection_Name = names(all_intersections_info),
    Count = unlist(all_intersections_info)
  )
  rownames(results_df) <- NULL # Clean up row names
  return(results_df)
}
