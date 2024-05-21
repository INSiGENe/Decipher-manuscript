#' Calculate Mean Decipher Scores by Cluster
#'
#' This function calculates the mean decipher score for each interaction within a cluster.
#' It also joins additional information (ligand and receptor) to the resulting data frame.
#'
#' @param scoresByRegulonCluster A data frame containing 'interaction', 'decipher_score',
#'        'ligand', and 'receptor' columns. It represents scores by regulon and cluster.
#'
#' @return A data frame with the mean decipher score for each interaction,
#'         sorted in descending order of scores. The returned data frame includes
#'         'interaction', 'decipher_score', 'ligand', 'receptor', and a new column
#'         'sender_cluster' with a fixed value 'mixed'.
#' @importFrom dplyr select group_by summarize ungroup left_join mutate arrange
#' @examples
#' @export
#' # Assuming 'scoresDF' is a data frame with the necessary columns:
#' meanScores <- calculateScoresByCluster(scoresDF)
calculateScoresByCluster <- function(scores_by_regulon_cluster){
  scores_by_regulon_cluster %>%
    select(interaction,decipher_score) %>%
    group_by(interaction) %>%
    summarize(decipher_score = mean(decipher_score)) %>%
    ungroup() %>%
    left_join(unique(select(scores_by_regulon_cluster,interaction,ligand,receptor)),by = "interaction")%>%
    mutate(sender_cluster = "mixed")%>%
    arrange(desc(decipher_score))
}
