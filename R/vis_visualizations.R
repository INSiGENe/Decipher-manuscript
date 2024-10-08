#' Pairwise Comparison of Prioritisation Scores
#'
#' This function takes in two data frames, conducts an inner join based on
#' interaction and receiver columns, calculates the Spearman correlation, and
#' produces a ggplot object displaying the correlation.
#'
#' @param df1 A data frame containing prioritisation scores from the first method.
#' @param df2 A data frame containing prioritisation scores from the second method.
#' @param method1 A character string representing the name of the first method.
#' @param method2 A character string representing the name of the second method.
#' @param overlap_max An integer representing a cap on the number of overlapping interactions such that the top overlap_max interactions from each method are selected
#'
#' @return A png saved in the figures folder displaying the correlation of prioritisation scores.
#' @export
#'
#' @importFrom graphics png dev.off
#' @importFrom stats cor rank
#' @importFrom utils file.path
#' @importFrom base dir.exists dir.create
#' @importFrom dplyr mutate arrange slice_head distinct inner_join bind_rows
#' @importFrom ggplot2 ggplot geom_point labs theme_minimal theme aes
#' @import ggplot2
#'
#' @examples
#' # pairwise_comparison_plot(nichenet_for_correlation, decipher_pre_processed, "NicheNet", "Decipher")
pairwise_comparison_plot <- function(df1, df2, method1, method2,overlap_max = NA,use_rank=FALSE,re_rank=FALSE){

  # Check if "figures" directory exists
  if (!dir.exists("figures")) {
    # Print a warning message
    warning("The 'figures' folder does not exist. Creating one now.")

    # Create the directory
    dir.create("figures")
  }

  if(use_rank){
    df1 <- df1 %>%
      mutate(prioritization_score = rank(-prioritization_score, ties.method = "min"))

    df2 <- df2 %>%
      mutate(prioritization_score = rank(-prioritization_score, ties.method = "min"))
  }


  # Inner join the data frames
  merged_data <- inner_join(df1, df2, by = c("interaction", "receiver"))

  if(!is.na(overlap_max)){


    library(dplyr)

    # For prioritization_score.x
    top_x <- merged_data %>%
      arrange(desc(prioritization_score.x)) %>%
      slice_head(n = overlap_max)

    # For prioritization_score.y
    top_y <- merged_data %>%
      arrange(desc(prioritization_score.y)) %>%
      slice_head(n = overlap_max)

    if(use_rank){
      # For prioritization_score.x
      top_x <- merged_data %>%
        arrange(prioritization_score.x) %>%
        slice_head(n = overlap_max)

      # For prioritization_score.y
      top_y <- merged_data %>%
        arrange(prioritization_score.y) %>%
        slice_head(n = overlap_max)
    }

    # Merging the two results
    merged_data <- bind_rows(top_x, top_y) %>%
      distinct() # Remove any potential duplicates

    if(use_rank & re_rank){
      merged_data <- merged_data %>%
        mutate(prioritization_score.x = rank(prioritization_score.x, ties.method = "min"),
               prioritization_score.y = rank(prioritization_score.y, ties.method = "min"))
    }

  }


  # Calculate the spearman correlation
  spearman_cor <- cor(merged_data$prioritization_score.x, merged_data$prioritization_score.y, method = "spearman")

  # Create the plot
  p <- ggplot(merged_data, aes(x = prioritization_score.x, y = prioritization_score.y)) +
    geom_point() +
    labs(
      title = paste("Prioritisation Score Correlation (", method1, " - ", method2, ")", sep = ""),
      subtitle = paste("Spearman Correlation: ", round(spearman_cor, 2)),
      x = method1,
      y = method2
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  filename <- paste(method1,method2,"correlation plot.png")
  filepath <- file.path("figures",filename)
  png(filepath,width=15,height=15,units="cm",res = 500)
  print(p)
  dev.off()
}

#' Plot Spearman Correlation Matrix
#'
#' Generates a heatmap from a Spearman correlation matrix for interaction data between methods.
#' The function modifies the input matrix to zero out the lower triangular part for better
#' visualization of the symmetric matrix and saves the heatmap as a PNG file.
#'
#' @param spearman_matrix A square matrix containing Spearman correlation coefficients
#' between different methods. It is assumed that the matrix is symmetric.
#' @param dataset_name A character string representing the name of the dataset,
#' which will be used as part of the filename for the saved heatmap image.
#' @param file_path The directory path where the heatmap image file will be saved.
#'
#' @details The function zeros out the lower triangular part of the spearman_matrix
#' to focus on the upper triangular part in the heatmap. It then saves the heatmap
#' as a PNG file named using the `dataset_name` in the specified `file_path`.
#' The heatmap visualizes the absolute values of the Spearman correlation coefficients
#' with color gradients.
#'
#' @return Invisibly returns NULL, as the main purpose is side-effect (image file creation).
#'
#' @examples
#' \dontrun{
#'   # Assuming `spearman_matrix` is your Spearman correlation matrix
#'   plotInteractionCorrelation(spearman_matrix, "MyDataset", "/path/to/save")
#' }
#'
#' @importFrom grDevices png dev.off
#' @importFrom gplots heatmap.2 colorpanel
#' @importFrom utils file.path
plotInteractionCorrelation <- function(spearman_matrix,dataset_name,file_path){
  spearman_matrix[lower.tri(spearman_matrix)] <- 0
  spearman_file <- paste(dataset_name,"spearman_matrix.png")
  png(file.path(file_path,spearman_file),width = 15,height = 15, units = "cm",res=600)
  p <- heatmap.2(abs(spearman_matrix),
                 trace="none",
                 dendrogram = "none",
                 Rowv = FALSE,
                 Colv = FALSE,
                 cexRow = 1.5,
                 cexCol = 1.5,
                 col = colorpanel(20, "white", "yellowgreen", "seagreen"),
                 key = TRUE,
                 keysize = 1.5,
                 key.xlab = "spear. cor.",
                 key.title = "",
                 key.ylab = "count",
                 margins = c(10,10))
  print(p)
  dev.off()
}

#' Plot Search Space Matrix
#'
#' Creates a heatmap visualization of the search space (k) matrix, which indicates
#' the minimum number of top-ranked interactions needed to achieve a specified number
#' of overlaps between pairs of methods. The lower triangle of the matrix is masked
#' to NA to avoid redundancy in the symmetric matrix. The heatmap is saved as a PNG
#' file.
#'
#' @param k_matrix A symmetric matrix where each element represents the minimum number
#' of top interactions required to achieve at least 100 overlapping interactions
#' between the methods corresponding to the row and column.
#' @param dataset_name A character string representing the name of the dataset, used
#' to generate the filename of the output PNG file.
#' @param file_path The directory path where the output PNG file will be saved.
#'
#' @return A heatmap visualization of the k matrix is saved as a PNG file at the
#' specified location. The function itself invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#'   plotSearchSpace(k_matrix, "my_dataset", "/path/to/save/directory")
#' }
#'
#' @importFrom grDevices png dev.off
#' @importFrom gplots heatmap.2 colorpanel
#' @importFrom utils file.path
plotSearchSpace <- function(k_matrix,dataset_name,file_path){
  k_matrix[lower.tri(k_matrix)] <- NA
  k_file <- paste(dataset_name,"k_matrix.png")
  png(file.path(file_path,k_file),width = 15,height = 15, units = "cm",res=600)
  p <- heatmap.2(k_matrix,
                 trace="none",
                 dendrogram = "none",
                 Rowv = FALSE,
                 Colv = FALSE,
                 cexRow = 1.5,
                 cexCol = 1.5,
                 col = colorpanel(20, "white", "tomato1", "tomato4"),
                 key = TRUE,
                 keysize = 1.5,
                 key.xlab = "k-value",
                 key.title = "",
                 key.ylab = "count",
                 margins = c(10,10))
  print(p)
  dev.off()
}

#' Plot UpSet Plot of Overlap Among Methods
#'
#' Generates an UpSet plot visualizing the overlap of top interactions among multiple
#' interaction prediction methods. The plot highlights the intersection and unique
#' contributions of each method to the pooled set of top interactions.
#'
#' @param method_results_lists A list of vectors, where each vector contains the top
#' interactions (as unique identifiers) predicted by a different method. The list should
#' be named with the names corresponding to the methods (e.g., 'NicheNet', 'LIANA+', etc.).
#' @param dataset_name A character string that specifies the dataset name, used to
#' generate the filename for the output plot.
#' @param file_path The path where the output plot PNG file should be saved.
#'
#' @return Saves an UpSet plot as a PNG file at the specified location. The function
#' itself invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#'   method_results <- list(
#'     NicheNet = c("interaction1", "interaction2"),
#'     `LIANA+` = c("interaction2", "interaction3"),
#'     Decipher = c("interaction1", "interaction3", "interaction4"),
#'     Connectome = c("interaction5"),
#'     NATMI = c("interaction2", "interaction4", "interaction5")
#'   )
#'   plotUpsetPlot(method_results, "my_dataset", "/path/to/save")
#' }
#'
#' @importFrom graphics png dev.off
#' @importFrom utils file.path paste
#' @importFrom UpSetR upset fromList intersects
plotUpsetPlot <- function(method_results_lists,dataset_name,file_path){
  upset_file <- paste(dataset_name,"upset_plot_overlap.png")
  png(file.path(file_path,upset_file),width=15,height=10,units="cm",res=600)
  p <- UpSetR::upset(UpSetR::fromList(method_results_lists),
                     sets = c("NicheNet", "LIANA+", "Decipher","Connectome","NATMI"),
                     sets.x.label = "Top 100 interactions",
                     mainbar.y.label = "Matching interactions",
                     text.scale = c(1.5,1.5,0.1,0.1,1.5,1.5),
                     queries = list(
                       list(query = UpSetR::intersects, params = list("Decipher"), color = "royalblue3", active = T)))

  print(p)
  dev.off()
}

#' Plot Cytosig Significance Matrix
#'
#' Visualizes a heatmap of the Cytosig significance matrix. This function processes
#' a dataframe containing Cytosig significance data, extracting unique values and
#' creating a matrix visualization. Significant values are highlighted with asterisks,
#' and the top ten rows with the highest absolute maximum values are displayed.
#' The heatmap is saved as a PNG file.
#'
#' @param cytosig_significance A dataframe containing Cytosig significance data,
#' expected to include columns for ligands, genes, and their respective significance
#' scores across various cell types or conditions.
#' @param output_figures_filepath A character string specifying the directory path
#' where the output PNG file should be saved.
#'
#' @return Invisibly returns NULL. A heatmap plot is generated and saved as a PNG file
#' at the specified location. The plot visualizes the top thirty rows based on the highest
#' absolute significance scores, with significance levels indicated by asterisks.
#'
#' @examples
#' \dontrun{
#'   plotCytosigSignificanceMatrix(cytosig_significance, "path/to/output/")
#' }
#'
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_fill_gradient2 theme_minimal theme labs scale_y_discrete
#' @importFrom dplyr select filter
#' @importFrom reshape2 melt
#' @importFrom graphics png dev.off
#' @importFrom utils file.path
plotCytosigSignificanceMatrix <- function(cytosig_significance,output_figures_filepath){
  cytosig_significance_matrix <- cytosig_significance %>% select(-ligand,-gene) %>% unique()
  cytosig_significance_matrix <- as.matrix(cytosig_significance_matrix)
  rownames(cytosig_significance_matrix) <- cytosig_significance$ligand %>% unique()

  # get the absolute maximum for each row
  abs_max_values <- apply(cytosig_significance_matrix, 1, get_abs_max)
  top_ten_rows <- names(abs_max_values)[order(abs_max_values,decreasing=TRUE)][1:30]

  color_breaks <- seq(-15, 15, length.out=20)

  # Convert the matrix to a long format data frame
  long_data <- reshape2::melt(cytosig_significance_matrix)

  # Add a column to determine the number of asterisks
  long_data$asterisks <- ifelse(abs(long_data$value) > 3, "**",
                                ifelse(abs(long_data$value) > 2, "*", ""))

  long_data <- long_data %>%
    filter(Var1 %in% top_ten_rows)

  #abs_max <- max(abs(cytosig_significance_matrix))
  abs_max <- 15

  png(file.path(output_figures_filepath,"cytosig_median_z.png"),width = 24,height=16,units="cm",res = 400)
  p <- ggplot(long_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = asterisks), color = "black", vjust = 0.5, fontface = "bold", size = 5) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limits = c(-abs_max, abs_max), space = "Lab",
                         name = "med. z score") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12),
      axis.text.y = element_text(angle = 0),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y.right = element_text(angle = 0,size=12),
      #axis.ticks.y.right = element_line(),
      legend.position = "bottom"
    ) +
    labs(x = "", y = '', title = '') +
    scale_y_discrete(position = "right")

  print(p)
  dev.off()

}

#' Custom Fill Color Scale for ggplot2
#'
#' Provides a custom gradient fill color scale for ggplot2 visualizations, transitioning
#' from dark orchid to white to dark green. This scale is particularly useful for
#' visualizations that require a diverging color scheme centered around a midpoint.
#'
#' @return A `scale_fill_gradient2` ggplot2 scale function, configured with a
#' custom set of colors and options.
#'
#' @examples
#' \dontrun{
#'   ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, fill = Sepal.Length)) +
#'     geom_tile() +
#'     scaleFillColor()
#' }
#'
#' @importFrom ggplot2 scale_fill_gradient2
#' @export
scaleFillColor <- function(){
  scale_fill_gradient2(
    low = "darkorchid4",
    mid = "white",
    high = "darkgreen",
    midpoint = 0,
    space = "Lab",
    na.value = "white",
    guide = "colourbar",
    aesthetics = "fill"
  )
}


#' Plot Spearman Correlation Heatmap
#'
#' Creates a heatmap visualization of the Spearman correlation matrix. If labels are
#' provided, they are used to annotate the heatmap; otherwise, a simpler heatmap is
#' generated. The heatmap is saved as a PNG file in a predetermined directory, with
#' the filename reflecting whether labels were used.
#'
#' @param spearman_matrix A matrix of Spearman correlation coefficients.
#' @param dataset A character string indicating the dataset name, used in generating
#' the output file name.
#' @param labels (Optional) A matrix of labels corresponding to the rows and columns
#' of the `spearman_matrix`. If provided, these labels are used to annotate the heatmap
#' cells. Note: This parameter is implicitly inferred from the context but not explicitly
#' defined in the function signature.
#'
#' @details The function generates two types of heatmaps based on the presence of
#' `labels`. If `labels` is `NULL`, a heatmap without cell annotations is created. If
#' `labels` are provided, they are assumed to come from a `label_matrix` (which must be
#' available in the scope where the function is called) and are used to annotate the
#' heatmap cells.
#'
#' The output PNG file is saved in the "figures" directory, with the filename incorporating
#' the dataset name and indicating whether the heatmap includes top overlapping
#' transcription factors (TFs).
#'
#' @examples
#' \dontrun{
#'   # Assuming 'spearman_mat' is your Spearman correlation matrix
#'   plotSpearmanHeatmap(spearman_mat, "my_dataset")
#'   # To plot with labels, ensure 'label_matrix' is defined in your environment
#'   plotSpearmanHeatmap(spearman_mat, "my_dataset", labels=TRUE)
#' }
#'
#' @importFrom grDevices png dev.off
#' @importFrom gplots heatmap.2 colorpanel
#' @importFrom utils file.path paste
plotSpearmanHeatmap <- function(spearman_matrix,labels=NULL,output_filepath,manuscript=FALSE){
  if(is.null(labels)){

    if(manuscript){
      filename <- "figure_3_panel_a1_spearman_regulon_scores.png"
    } else {
      filename <- paste("spearman_regulon_scores.png")
    }

    png(file.path(output_filepath,filename),width = 15,height = 15, units = "cm",res=400)
    p <- heatmap.2(abs(spearman_matrix),
                   trace="none",
                   dendrogram = "none",
                   Rowv = FALSE,
                   Colv = FALSE,
                   cexRow = 0.7,
                   cexCol = 0.7,
                   col = colorpanel(20, "white", "yellowgreen", "seagreen"),
                   key = TRUE,
                   keysize = 1.5,
                   key.xlab = "spear. cor.",
                   key.title = "",
                   key.ylab = "count",
                   margins = c(8,8),
                   xlab = "CT",
                   ylab = "CT")
    print(p)
    dev.off()
  }else{
    if(manuscript){
      filename <- "figure_3_panel_a1_spearman_regulon_scores_top_overlapping_tf.png"
    } else {
      filename <- paste("spearman_regulon_scores_top_overlapping_tf.png")
    }
    png(file.path(output_filepath,filename),width = 15,height = 15, units = "cm",res=400)
    p <- heatmap.2(abs(spearman_matrix),
                   trace="none",
                   dendrogram = "none",
                   Rowv = FALSE,
                   Colv = FALSE,
                   cexRow = 0.7,
                   cexCol = 0.7,
                   col = colorpanel(20,"white", "yellowgreen", "seagreen"),
                   key = TRUE,
                   keysize = 1.5,
                   key.xlab = "spear. cor.",
                   key.title = "",
                   key.ylab = "count",
                   margins = c(8,8),
                   xlab = "CT",
                   ylab = "CT",
                   cellnote = label_matrix[rownames(spearman_matrix),colnames(spearman_matrix)],
                   notecol="black", notecex=0.9)
    print(p)
    dev.off()
  }
}



#' Create a bubble plot with gradient fill and customizable aesthetics
#'
#' This function generates a bubble plot using ggplot2, with various customizable aesthetic features.
#' It plots an interaction variable against a specified x-axis variable, coloring the points based on a color variable and scaling by size.
#' Additional customizations include gradient color filling, position of the plot, and optional axis adjustments.
#'
#' @param df Data frame containing the variables for plotting.
#' @param color.var A string specifying the column in `df` used for color gradient.
#' @param size.var A string specifying the column in `df` used for point size.
#' @param stroke.var A string specifying the stroke width of points.
#' @param plot.position A string indicating the plot position ('middle', 'right', or other).
#' @param col.min.val The minimum value for the color gradient scale.
#' @param col.max.val The maximum value for the color gradient scale.
#' @param plot.title String for the plot's title.
#' @param x_lab Label for the x-axis.
#' @param y_lab Label for the y-axis.
#' @param x_var A string specifying the x-axis variable.
#' @return A ggplot object representing the bubble plot.
#' @examples
#' df <- data.frame(interaction = rnorm(10), A = rnorm(10), B = sample(1:100, 10), C = sample(1:5, 10, replace = TRUE))
#' plotBubble(df, color.var="A", size.var="B", stroke.var="C", plot.position="middle",
#'            col.min.val=-2, col.max.val=2, plot.title="My Bubble Plot", x_lab="X Axis", y_lab="Y Axis", x_var="A")
#' @import ggplot2
#' @export
plotBubble <- function(df,color.var,size.var,stroke.var,plot.position,col.min.val,col.max.val,plot.title,x_lab,y_lab,x_var){

  this.plot <- ggplot(df,aes_string(y="interaction", x=x_var,fill = color.var)) +
    geom_point(aes_string(size = size.var,stroke=stroke.var), shape = 21) +
    labs(x = x_lab, y = y_lab) +
    theme_bw()+
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(col.min.val,col.max.val)
    )+ggtitle(label = plot.title)+ guides(size = "none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size=5),
          axis.text.y = element_text(size=5),
          legend.position = "bottom",
          legend.key.size = unit(0.4, 'cm'),           # Smaller legend keys
          legend.text = element_text(size = 6),        # Smaller legend text
          plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 2, unit = "pt"),
          plot.title = element_text(hjust = 0.5))

  if(plot.position == "middle"){
    this.plot <- this.plot + theme(axis.text.y=element_blank(),
                                   axis.ticks.y=element_blank())
  } else if (plot.position == "right"){
    this.plot <- this.plot + scale_y_discrete(position = "right")+
      theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 10, unit = "pt"))
  }
  return(this.plot)
}


#' Plot Decipher Prioritized Map
#'
#' This function reads and processes data to create and save a prioritized interaction map
#' for ligand-receptor interactions using the Decipher package.
#'
#' @param dataset_path Character. The path to the dataset directory containing the data files.
#'
#' @return This function does not return a value. It saves a plot and a CSV file in the specified dataset path.
#'
#' @examples
#' \dontrun{
#' plotDecipherPrioritizedMap("path/to/dataset")
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom grDevices png dev.off
#' @importFrom utils write.csv
#' @import patchwork
plotDecipherPrioritizedMap <- function(dataset_path,top_n,selected_receiver_cells = NULL,sc_feature_statistics=FALSE,primary_ct = NULL,split_by_direction = FALSE){
  decipher_path <- file.path(dataset_path,"data")
  #read data ----
  lr_markers_by_cluster <- readRDS(file.path(decipher_path,"lr_markers_by_cluster.rds"))
  feature_statistics <- readRDS(file.path(decipher_path,"feature_statistics.rds"))
  decipher_scores_by_cluster <- readRDS(file.path(decipher_path,"decipher_scores_by_cluster.rds"))
  L_set <- readRDS(file.path(decipher_path,"L_set.rds"))
  if(sc_feature_statistics){
    l_markers_by_cluster_sc <- readRDS(file.path(decipher_path,"l_markers_by_cluster_sc.rds"))
    l_markers_by_cluster_sc <- getLigandReceptorDiffExprMarkersByCt(l_markers_by_cluster_sc)

  }

  ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)
  decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)

  #normalize feature statistics
  normalized_feature_statistics <- feature_statistics
  normalized_feature_statistics$normalized.counts <- normalized_feature_statistics$sum.counts/normalized_feature_statistics$n.cell
  #total counts feature condition
  normalized_feature_statistics <- normalized_feature_statistics %>%
    group_by(condition,feature) %>%
    mutate(total.normalized.counts = sum(normalized.counts)) %>%
    ungroup() %>%
    mutate(frac.normalized.counts.features.condition = normalized.counts/total.normalized.counts)


  #PANEL D ----
  ##data wrangling ----
  #first we enrich decipher results with the information we will need for downstream visualization
  ##enrich ----
  #used to be called decipher_scores_by_cluster_df_enriched
  decipher_scores_by_cluster_bound_enriched <- decipher_scores_by_cluster_bound %>%
    select(interaction,receiver_cluster,decipher_score) %>%
    # Ensure complete combinations of interaction and receiver_cluster.
    tidyr::complete(interaction,receiver_cluster) %>%
    # Replace missing decipher scores with 0 and add sender_cluster column with value "mixed".
    mutate(decipher_score =tidyr::replace_na(decipher_score,0),
           sender_cluster = "mixed") %>%
    # Left join with ligand-receptor interaction data.
    left_join(select(L_set,ligand,receptor,interaction),by = "interaction") %>%
    # Add differential expression data for ligands.
    left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("receiver_cluster"="cluster","ligand"="gene")) %>%
    rename(ligand.diff.expr = avg_log2FC) %>%
    # Add differential expression data for receptors.
    left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("receiver_cluster"="cluster","receptor"="gene")) %>%
    rename(receptor.diff.expr = avg_log2FC) %>%
    # Replace missing values in differential expression columns with 0 and add condition column with value "case".
    mutate(ligand.diff.expr = replaceNAw0(ligand.diff.expr),
           receptor.diff.expr = replaceNAw0(receptor.diff.expr),
           condition = "case") %>%
    # Left join with normalized feature statistics for ligands.
    left_join(select(normalized_feature_statistics,cluster,feature,condition,frac.normalized.counts.features.condition), by = c("receiver_cluster"="cluster","ligand"="feature","condition"))%>%
    rename(ligand.frac = frac.normalized.counts.features.condition) %>%
    # Left join with normalized feature statistics for receptors.
    left_join(select(normalized_feature_statistics,cluster,feature,condition,frac.normalized.counts.features.condition), by = c("receiver_cluster"="cluster","receptor"="feature","condition"))%>%
    rename(receptor.frac = frac.normalized.counts.features.condition)

  ## top interactions ----
  #used to be called decipher_for_overlap/merged_data and was comprised of decipher_bound and decipher_pre_processed
  decipher_scores_by_cluster_bound_clean <- decipher_scores_by_cluster_bound %>%
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    rename(receiver=receiver_cluster,sender=sender_cluster,prioritization_score=decipher_score) %>%
    select(interaction,ligand,receptor,receiver,prioritization_score) %>%
    arrange(prioritization_score)


  #Merge the results from the three methods
  # selected_cts <- c("B_cell","Monocyte","CD4_T","NK_cell_1")
  # selected_rcts <- c("B_cell","Monocyte","CD4_T","NK_cell_1")
  # selected_scts <- c("B_cell","Monocyte","HSC","CD4_T","CD8_Tem","NK_cell_1")

  #now
  if(is.null(selected_receiver_cells)){
    if(!split_by_direction){
      decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
        mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
        group_by(decipher_score_sign) %>%
        arrange(desc(abs(prioritization_score))) %>%
        select(interaction) %>%
        distinct() %>%
        slice_head(n = top_n) %>%
        ungroup() %>%
        left_join(decipher_scores_by_cluster_bound)
    } else {
      decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
        mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
        filter(decipher_score_sign == "positive") %>%
        arrange(desc(abs(prioritization_score))) %>%
        select(interaction) %>%
        distinct() %>%
        slice_head(n = top_n) %>%
        left_join(decipher_scores_by_cluster_bound)
    }

  } else {
    if(is.null(primary_ct)){
      if(!split_by_direction){
        decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
          mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
          filter(receiver %in% selected_receiver_cells) %>%
          group_by(decipher_score_sign) %>%
          arrange(desc(abs(prioritization_score))) %>%
          select(interaction) %>%
          distinct() %>%
          slice_head(n = top_n) %>%
          ungroup() %>%
          left_join(decipher_scores_by_cluster_bound)
      } else {
        decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
          mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
          filter(receiver %in% selected_receiver_cells) %>%
          filter(decipher_score_sign == "positive") %>%
          arrange(desc(abs(prioritization_score))) %>%
          select(interaction) %>%
          distinct() %>%
          slice_head(n = top_n) %>%
          left_join(decipher_scores_by_cluster_bound)
      }
    } else {
        if(!split_by_direction){
          decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
            mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
            filter(receiver %in% primary_ct) %>%
            group_by(decipher_score_sign) %>%
            arrange(desc(abs(prioritization_score))) %>%
            select(interaction) %>%
            distinct() %>%
            slice_head(n = top_n) %>%
            ungroup() %>%
            left_join(decipher_scores_by_cluster_bound)
        } else {
          decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
            mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
            filter(receiver %in% primary_ct) %>%
            filter(decipher_score_sign == "positive") %>%
            arrange(desc(abs(prioritization_score))) %>%
            select(interaction) %>%
            distinct() %>%
            slice_head(n = top_n) %>%
            left_join(decipher_scores_by_cluster_bound)
        }

    }


  }



  top_interactions <- decipher_top_interactions_all_clusters %>%
    select(interaction) %>%
    distinct() %>%
    unlist(use.names=FALSE)

  ## visualization ----
  base_data <- decipher_scores_by_cluster_bound_enriched %>%
    filter(interaction %in% top_interactions) %>%
    mutate(size = 1)

  plot_limits_ligand <- list(max = max(base_data$ligand.diff.expr),min = min(base_data$ligand.diff.expr))

  base_data$stroke <- 0.5

  base_data <- base_data %>%
    mutate(stroke_ligand = if_else(ligand.frac > 0.05,0.5,NA)) %>%
    mutate(size_ligand = if_else(ligand.frac > 0.05,ligand.frac,NA)) %>%
    mutate(stroke_receptor = if_else(receptor.frac > 0.05,0.5,NA),
           size_receptor = if_else(receptor.frac > 0.05,receptor.frac,NA)) %>%
    mutate(size = if_else(abs(decipher_score) > 0.1,1,NA))%>%
    mutate(receiver_cluster=if_else(receiver_cluster == "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells","C8",receiver_cluster))


  base_data$receiver_cluster <- convert_text_patterns(base_data$receiver_cluster)

  if(!sc_feature_statistics){
    base_data_ligand <- base_data
  } else {
    base_data_ligand <- normalized_feature_statistics %>%
      filter(condition == "case") %>%  # Assuming 'case' condition is desired
      select(cluster, feature, frac.normalized.counts.features.condition) %>%
      filter(feature %in% unique(base_data$ligand)) %>%
      # Join with base_data to get additional information
      rename(ligand.frac = frac.normalized.counts.features.condition) %>%
      # Handle NAs by setting them to 0 or appropriate default
      mutate(
        ligand.frac = replace_na(ligand.frac, 0),
      ) %>%
      left_join(l_markers_by_cluster_sc[,c("cluster","gene","avg_log2FC")],by = c("cluster"="cluster","feature"="gene")) %>%
      rename(ligand.diff.expr = avg_log2FC) %>%
      left_join(unique(select(base_data,interaction,ligand)),by = c("feature"="ligand"))%>%
      dplyr::rename("receiver_cluster"="cluster") %>%
      dplyr::rename("ligand"="feature")

    plot_limits_ligand <- list(max = max(base_data_ligand$ligand.diff.expr),min = min(base_data_ligand$ligand.diff.expr))

    base_data_ligand$stroke <- 0.5

    base_data_ligand <- base_data_ligand %>%
      mutate(stroke_ligand = if_else(ligand.frac > 0.05,0.5,NA)) %>%
      mutate(size_ligand = if_else(ligand.frac > 0.05,ligand.frac,NA))
  }


  ligand_bubble_plot <- plotBubble(
    df = base_data_ligand,
    x_var = "receiver_cluster" ,
    color.var = "ligand.diff.expr",
    size.var = "size_ligand",
    stroke.var = "stroke_ligand",
    plot.position = "left",
    col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
    plot.title = "Ligand",
    x_lab= "SCT",
    y_lab = "Interaction")

  if(is.null(selected_receiver_cells)){
    base_data_decipher <- base_data

  } else {
    base_data_decipher <- base_data %>%
      filter(receiver_cluster %in% selected_receiver_cells)
  }

  plot_limits_receptor <- list(max = max(base_data_decipher$receptor.diff.expr),min = min(base_data_decipher$receptor.diff.expr))
  plot_limits_decipher <- list(max = max(base_data_decipher$decipher_score),min = min(base_data_decipher$decipher_score))

  decipher_bubble_plot <- plotBubble(
    df = base_data_decipher,
    x_var = "receiver_cluster",
    color.var = "decipher_score",
    size.var = "size",
    stroke.var = "stroke",
    plot.position = "middle",
    col.min.val=plot_limits_decipher$min,col.max.val=plot_limits_decipher$max,
    plot.title = "Decipher score",
    x_lab= "RCT",
    y_lab = "")

  receptor_bubble_plot <- plotBubble(
    df = base_data_decipher,
    x_var = "receiver_cluster",
    color.var = "receptor.diff.expr",
    size.var = "size_receptor",
    stroke.var = "stroke_receptor",
    plot.position = "middle",
    col.min.val=plot_limits_receptor$min,col.max.val=plot_limits_receptor$max,
    plot.title = "Receptor",
    x_lab= "RCT",
    y_lab = "")

  if(is.null(primary_ct)){
    if(!split_by_direction){
      png(file.path(dataset_path,"figures","decipher_plot_prioritized.png"),width = 21,height = 11,units = "cm",res = 600)

    } else {
      png(file.path(dataset_path,"figures","decipher_plot_prioritized_split.png"),width = 21,height = 11,units = "cm",res = 600)

    }
    print(ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot+
            patchwork::plot_layout(widths = c(2, 1, 1)))
    dev.off()

    write.csv(base_data,file.path(dataset_path,"figures","decipher_plot_prioritized.csv"))

  } else {
    if(!split_by_direction){
      filename_no_ext <- paste("decipher_plot_prioritized_",primary_ct,sep="")
    } else {
      filename_no_ext <- paste("decipher_plot_prioritized_","split_",primary_ct,sep="")

    }
    filename <- paste(filename_no_ext,".png",sep="")
    png(file.path(dataset_path,"figures",filename),width = 21,height = 11,units = "cm",res = 600)
    print(ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot+
            patchwork::plot_layout(widths = c(2, 1, 1)))
    dev.off()

    write.csv(base_data,file.path(dataset_path,"figures",paste(filename_no_ext,".csv",sep="")))

  }

}

#' Plot Decipher Prioritized Map
#'
#' This function reads and processes data to create and save a prioritized interaction map
#' for ligand-receptor interactions using the Decipher package.
#'
#' @param dataset_path Character. The path to the dataset directory containing the data files.
#' @param receiver_cell_type Character. The receiver cell type to filter interactions by.
#' @param output_filename Character. The filename for the saved plot and CSV file.
#' @param log_transform log transform
#' @param slice_n number of top interactions
#' @param return_plot_object
#'
#' @return This function does not return a value. It saves a plot and a CSV file in the specified dataset path.
#'
#' @examples
#' \dontrun{
#' plotDecipherPrioritizedMap("path/to/dataset", "CD4_T", "decipher_plot_prioritized")
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom grDevices png dev.off
#' @importFrom utils write.csv
plotDecipherPrioritizedMap_v2 <- function(dataset_path, receiver_cell_type = NULL, output_filename,log_transform,slice_n = 5,return_plot_object = FALSE){
  decipher_path <- file.path(dataset_path, "data")

  # Read data ----
  lr_markers_by_cluster <- readRDS(file.path(decipher_path, "lr_markers_by_cluster.rds"))
  feature_statistics <- readRDS(file.path(decipher_path, "feature_statistics.rds"))
  decipher_scores_by_cluster <- readRDS(file.path(decipher_path, "decipher_scores_by_cluster.rds"))
  L_set <- readRDS(file.path(decipher_path, "L_set.rds"))

  ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)
  decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)

  # Normalize feature statistics
  normalized_feature_statistics <- feature_statistics
  normalized_feature_statistics$normalized.counts <- normalized_feature_statistics$sum.counts / normalized_feature_statistics$n.cell
  normalized_feature_statistics <- normalized_feature_statistics %>%
    group_by(condition, feature) %>%
    mutate(total.normalized.counts = sum(normalized.counts)) %>%
    ungroup() %>%
    mutate(frac.normalized.counts.features.condition = normalized.counts / total.normalized.counts)

  # Enrich decipher results ----
  decipher_scores_by_cluster_bound_enriched <- decipher_scores_by_cluster_bound %>%
    select(interaction, receiver_cluster, decipher_score) %>%
    tidyr::complete(interaction, receiver_cluster) %>%
    mutate(decipher_score = tidyr::replace_na(decipher_score, 0),
           sender_cluster = "mixed") %>%
    left_join(select(L_set, ligand, receptor, interaction), by = "interaction") %>%
    left_join(ct_lr_markers[, c("cluster", "gene", "avg_log2FC")], by = c("receiver_cluster" = "cluster", "ligand" = "gene")) %>%
    rename(ligand.diff.expr = avg_log2FC) %>%
    left_join(ct_lr_markers[, c("cluster", "gene", "avg_log2FC")], by = c("receiver_cluster" = "cluster", "receptor" = "gene")) %>%
    rename(receptor.diff.expr = avg_log2FC) %>%
    mutate(ligand.diff.expr = replaceNAw0(ligand.diff.expr),
           receptor.diff.expr = replaceNAw0(receptor.diff.expr),
           condition = "case") %>%
    left_join(select(normalized_feature_statistics, cluster, feature, condition, frac.normalized.counts.features.condition), by = c("receiver_cluster" = "cluster", "ligand" = "feature", "condition")) %>%
    rename(ligand.frac = frac.normalized.counts.features.condition) %>%
    left_join(select(normalized_feature_statistics, cluster, feature, condition, frac.normalized.counts.features.condition), by = c("receiver_cluster" = "cluster", "receptor" = "feature", "condition")) %>%
    rename(receptor.frac = frac.normalized.counts.features.condition)

  # Top interactions ----
  decipher_scores_by_cluster_bound_clean <- decipher_scores_by_cluster_bound %>%
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    rename(receiver = receiver_cluster, sender = sender_cluster, prioritization_score = decipher_score) %>%
    select(interaction, ligand, receptor, receiver, prioritization_score) %>%
    arrange(prioritization_score) %>%
    mutate(decipher_score_sign = if_else(prioritization_score >= 0, "positive", "negative"))

if(is.null(receiver_cell_type)){
  decipher_top_interactions_by_receiver <- decipher_scores_by_cluster_bound_clean %>%
    #mutate(decipher_score_sign = if_else(prioritization_score >= 0, "positive", "negative")) %>%
    group_by(decipher_score_sign) %>%
    arrange(desc(abs(prioritization_score))) %>%
    select(interaction,decipher_score_sign) %>%
    distinct() %>%
    slice_head(n = slice_n) %>%
    ungroup() %>%
    left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction"))
} else {
  decipher_top_interactions_by_receiver <- decipher_scores_by_cluster_bound_clean %>%
    filter(receiver == receiver_cell_type) %>%
    #mutate(decipher_score_sign = if_else(prioritization_score >= 0, "positive", "negative")) %>%
    group_by(decipher_score_sign) %>%
    arrange(desc(abs(prioritization_score))) %>%
    select(interaction,decipher_score_sign) %>%
    distinct() %>%
    slice_head(n = slice_n) %>%
    ungroup() %>%
    left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction"))
}



  top_interactions <- decipher_top_interactions_by_receiver %>%
    select(interaction) %>%
    distinct() %>%
    unlist(use.names = FALSE)

  # Log-transform decipher scores
  if(log_transform){
    decipher_scores_by_cluster_bound_enriched <- decipher_scores_by_cluster_bound_enriched %>%
      mutate(decipher_score = case_when(
        decipher_score < -1e-10 ~ -log(abs(decipher_score)),
        decipher_score > 1e-10 ~ log(decipher_score),
        TRUE ~ 0
      ))
  }


  # Visualization ----
  base_data <- decipher_scores_by_cluster_bound_enriched %>%
    filter(interaction %in% top_interactions) %>%
    mutate(size = 1)

  plot_limits_ligand <- list(max = max(base_data$ligand.diff.expr), min = min(base_data$ligand.diff.expr))
  plot_limits_receptor <- list(max = max(base_data$receptor.diff.expr), min = min(base_data$receptor.diff.expr))
  plot_limits_decipher <- list(max = max(base_data$decipher_score), min = min(base_data$decipher_score))

  base_data$stroke <- 0.5

  base_data <- base_data %>%
    mutate(stroke_ligand = if_else(ligand.frac > 0.05, 0.5, NA)) %>%
    mutate(size_ligand = if_else(ligand.frac > 0.05, ligand.frac, NA)) %>%
    mutate(stroke_receptor = if_else(receptor.frac > 0.05, 0.5, NA),
           size_receptor = if_else(receptor.frac > 0.05, receptor.frac, NA)) %>%
    mutate(size = if_else(abs(decipher_score) > 0.1, 1, NA)) %>%
    mutate(receiver_cluster = if_else(receiver_cluster == "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells", "C8", receiver_cluster))

  base_data$receiver_cluster <- convert_text_patterns(base_data$receiver_cluster)

  ligand_bubble_plot <- plotBubble(
    df = base_data,
    x_var = "receiver_cluster",
    color.var = "ligand.diff.expr",
    size.var = "size_ligand",
    stroke.var = "stroke_ligand",
    plot.position = "left",
    col.min.val = plot_limits_ligand$min, col.max.val = plot_limits_ligand$max,
    plot.title = "Ligand",
    x_lab = "SCT",
    y_lab = "Interaction")

  decipher_bubble_plot <- plotBubble(
    df = base_data,
    x_var = "receiver_cluster",
    color.var = "decipher_score",
    size.var = "size",
    stroke.var = "stroke",
    plot.position = "middle",
    col.min.val = plot_limits_decipher$min, col.max.val = plot_limits_decipher$max,
    plot.title = "Decipher Score",
    x_lab = "RCT",
    y_lab = "")

  receptor_bubble_plot <- plotBubble(
    df = base_data,
    x_var = "receiver_cluster",
    color.var = "receptor.diff.expr",
    size.var = "size_receptor",
    stroke.var = "stroke_receptor",
    plot.position = "middle",
    col.min.val = plot_limits_receptor$min, col.max.val = plot_limits_receptor$max,
    plot.title = "Receptor",
    x_lab = "RCT",
    y_lab = "")

  output_plot_path <- file.path(dataset_path, "figures", paste0(output_filename, ".png"))
  output_csv_path <- file.path(dataset_path, "figures", paste0(output_filename, ".csv"))

  if(return_plot_object){
    p <-  ligand_bubble_plot + decipher_bubble_plot + receptor_bubble_plot
    return(p)
  }else{
    png(output_plot_path, width = 21, height = 11, units = "cm", res = 600)
    print(ligand_bubble_plot + decipher_bubble_plot + receptor_bubble_plot)
    dev.off()
    write.csv(base_data, output_csv_path)

  }


}


