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
      dplyr::slice_head(n = overlap_max)

    # For prioritization_score.y
    top_y <- merged_data %>%
      arrange(desc(prioritization_score.y)) %>%
      dplyr::slice_head(n = overlap_max)

    if(use_rank){
      # For prioritization_score.x
      top_x <- merged_data %>%
        arrange(prioritization_score.x) %>%
        dplyr::slice_head(n = overlap_max)

      # For prioritization_score.y
      top_y <- merged_data %>%
        arrange(prioritization_score.y) %>%
        dplyr::slice_head(n = overlap_max)
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
    labs(x = NULL, y = y_lab) +
    theme_bw()+
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
      aesthetics = "fill",
      limits = c(col.min.val,col.max.val)
    )+ggtitle(label = plot.title)+ guides(size = "none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, 'cm'),           # Smaller legend keys
          legend.text = element_text(size = 10),        # Smaller legend text
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
plotDecipherPrioritizedMap <- function(dataset_path,top_n,selected_receiver_cells = NULL,sc_feature_statistics=FALSE,primary_ct = NULL,split_by_direction = FALSE,direction = c("pos","neg"),dataset_name,abs_decipher_plot_limit = NULL,priority_receiver_cells = NULL,width=NULL,height=NULL
){
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
  decipher_top_interactions_all_clusters <- decipher_scores_by_cluster_bound_clean %>%
    mutate(decipher_score_sign = if_else(prioritization_score >= 0, "positive", "negative"))

  # Apply filter based on selected receiver cells or primary_ct
  if (!is.null(selected_receiver_cells)) {
    decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
      filter(receiver %in% selected_receiver_cells)
  }

  if (!is.null(primary_ct)) {
    decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
      filter(receiver %in% primary_ct)
  }

  # Further filtering based on direction if required
  if (split_by_direction) {
    if(direction == "pos"){
      decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
        filter(decipher_score_sign == "positive")
    } else {
      decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
        filter(decipher_score_sign == "negative")
    }

  }

  if (!is.null(priority_receiver_cells)) {
  decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
    mutate(priority = if_else(receiver %in% priority_receiver_cells, 1, 0)) %>%
    arrange(desc(priority), desc(abs(prioritization_score)))
} else {
  decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
    arrange(desc(abs(prioritization_score)))
}

    
    # Then continue normally
decipher_top_interactions_all_clusters <- decipher_top_interactions_all_clusters %>%
  group_by(decipher_score_sign) %>%
  select(interaction) %>%
  distinct() %>%
  dplyr::slice_head(n = top_n) %>%
  ungroup() %>%
  left_join(decipher_scores_by_cluster_bound)



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
  if(!is.null(selected_receiver_cells)){
    base_data$receiver_cluster <- factor(base_data$receiver_cluster,levels = selected_receiver_cells)

  }

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

  if (dataset_name == "sample_1") {
    base_data_ligand <- base_data_ligand %>%
      mutate(receiver_cluster = case_when(
        grepl("B", receiver_cluster) ~ "B",
        grepl("CD14", receiver_cluster) ~ "CD14+ M",
        grepl("CD4", receiver_cluster) ~ "CD4 T",
        grepl("CD8", receiver_cluster) ~ "CD8 T",
        TRUE ~ receiver_cluster  # Keep original if no match
      ))
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
  if(is.null(abs_decipher_plot_limit)){
    plot_limits_decipher <- list(max = max(base_data_decipher$decipher_score),min = min(base_data_decipher$decipher_score))
  } else{
    # Condition 2: Plot limit is specified
    plot_limits_decipher <- list(
      max = abs_decipher_plot_limit,
      min = -abs_decipher_plot_limit
    )    # Calculate epsilon

    epsilon <- 0.01 * abs_decipher_plot_limit

    # Update 'decipher_score' by capping the values
    base_data_decipher <- base_data_decipher %>%
      mutate(decipher_score = case_when(
        decipher_score > abs_decipher_plot_limit ~ abs_decipher_plot_limit - epsilon,
        decipher_score < -abs_decipher_plot_limit ~ -abs_decipher_plot_limit + epsilon,
        TRUE ~ decipher_score
      ))

  }

  if (dataset_name == "sample_1") {
    base_data_decipher <- base_data_decipher %>%
      mutate(receiver_cluster = case_when(
        grepl("B", receiver_cluster) ~ "B",
        grepl("CD14", receiver_cluster) ~ "CD14+ M",
        grepl("CD4", receiver_cluster) ~ "CD4 T",
        grepl("CD8", receiver_cluster) ~ "CD8 T",
        TRUE ~ receiver_cluster  # Keep original if no match
      ))
  }

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
  # Helper function to handle file saving
  save_plot_and_data <- function(filename_prefix, ligand_plot, decipher_plot, receptor_plot, base_data, dataset_path,width = NULL,height = NULL) {
    # Construct file paths
    filename_png <- paste0(filename_prefix, ".png")
    filename_csv <- paste0(filename_prefix, ".csv")

    # Save plot as PNG
    if(is.null(width) | is.null(height)){
          png(file.path(dataset_path, "figures", filename_png), width = 21, height = 11, units = "cm", res = 600)
    } else {
      png(file.path(dataset_path, "figures", filename_png), width = width, height = height, units = "cm", res = 600)
    }
    print(ligand_plot + decipher_plot + receptor_plot + patchwork::plot_layout(widths = c(2, 1, 1)))
    dev.off()

    # Save base data as CSV
    write.csv(base_data, file.path(dataset_path, "figures", filename_csv))
  }

  # Logic to determine the filename
  if (is.null(primary_ct)) {
    if (!split_by_direction) {
      filename_prefix <- paste(dataset_name,"decipher_plot_prioritized",sep="_")
    } else {
      if (direction == "pos") {
        filename_prefix <- paste(dataset_name,"decipher_plot_prioritized_split_pos",sep="_")
      } else {
        filename_prefix <- paste(dataset_name,"decipher_plot_prioritized_split_neg",sep="_")
      }
    }
  } else {
    if (!split_by_direction) {
      filename_prefix <- paste(dataset_name,"decipher_plot_prioritized", primary_ct,sep = "_")
    } else {
      if (direction == "pos") {
        filename_prefix <- paste(dataset_name,"decipher_plot_prioritized_split_pos", primary_ct,sep="_")
      } else {
        filename_prefix <- paste(dataset_name,"decipher_plot_prioritized_split_neg", primary_ct,sep="_")
      }
    }
  }

  # Call the helper function to save plots and data
  save_plot_and_data(filename_prefix, ligand_bubble_plot, decipher_bubble_plot, receptor_bubble_plot, base_data, dataset_path,width,height)


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
    dplyr::slice_head(n = slice_n) %>%
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
    dplyr::slice_head(n = slice_n) %>%
    ungroup() %>%
    left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction","receiver"))
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


plotDecipherPrioritizedMap_v3 <- function(dataset_path, receiver_cell_type = NULL, output_filename, log_transform, slice_n = 5, return_plot_object = FALSE) {
  library(ggplot2)
  library(patchwork)  # For improved layout management
  
  decipher_path <- file.path(dataset_path, "data")

  # Read data ----
  lr_markers_by_cluster <- readRDS(file.path(decipher_path, "lr_markers_by_cluster.rds"))
  feature_statistics <- readRDS(file.path(decipher_path, "feature_statistics.rds"))
  decipher_scores_by_cluster <- readRDS(file.path(decipher_path, "decipher_scores_by_cluster.rds"))
  L_set <- readRDS(file.path(decipher_path, "L_set.rds"))

  ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)
  decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)

  # Normalize feature statistics ----
  normalized_feature_statistics <- feature_statistics %>%
    mutate(normalized.counts = sum.counts / n.cell) %>%
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

  if (is.null(receiver_cell_type)) {
    decipher_top_interactions_by_receiver <- decipher_scores_by_cluster_bound_clean %>%
      group_by(decipher_score_sign) %>%
      arrange(desc(abs(prioritization_score))) %>%
      select(interaction, decipher_score_sign) %>%
      distinct() %>%
      dplyr::slice_head(n = slice_n) %>%
      ungroup() %>%
      left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction"))
  } else {
    decipher_top_interactions_by_receiver <- decipher_scores_by_cluster_bound_clean %>%
      filter(receiver == receiver_cell_type) %>%
      group_by(decipher_score_sign) %>%
      arrange(desc(abs(prioritization_score))) %>%
      select(interaction, decipher_score_sign) %>%
      distinct() %>%
      dplyr::slice_head(n = slice_n) %>%
      ungroup() %>%
      left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction"))
  }

  top_interactions <- decipher_top_interactions_by_receiver %>%
    select(interaction) %>%
    distinct() %>%
    unlist(use.names = FALSE)

  # Log-transform decipher scores
  if (log_transform) {
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

  # Define color scales and limits
  color_scale <- scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

  base_plot <- ggplot(base_data, aes(x = receiver_cluster, y = interaction)) +
    geom_point(aes(size = abs(decipher_score), fill = decipher_score), shape = 21, color = "black", stroke = 0.3) +
    color_scale +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom") +
    labs(x = "Receiver Cluster", y = "", title = "Decipher Prioritized Interactions")+
    guides(
        fill = guide_colorbar(order = 1),  # Gradient legend on top
        size = guide_legend(order = 2)     # Size legend below it
    )

  # Ligand Expression Plot
  ligand_plot <- ggplot(base_data, aes(x = receiver_cluster, y = interaction)) +
    geom_point(aes(size = abs(ligand.diff.expr), fill = ligand.diff.expr), shape = 21, color = "black", stroke = 0.3) +
    color_scale +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom") +
    labs(x = "Receiver Cluster", y = "Interaction", title = "Ligand Expression")+
    guides(
        fill = guide_colorbar(order = 1),  # Gradient legend on top
        size = guide_legend(order = 2)     # Size legend below it
    )

  # Receptor Expression Plot
  receptor_plot <- ggplot(base_data, aes(x = receiver_cluster, y = interaction)) +
    geom_point(aes(size = abs(receptor.diff.expr), fill = receptor.diff.expr), shape = 21, color = "black", stroke = 0.3) +
    color_scale +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom") +
    labs(x = "Receiver Cluster", y = "", title = "Receptor Expression")+
    guides(
        fill = guide_colorbar(order = 1),  # Gradient legend on top
        size = guide_legend(order = 2)     # Size legend below it
    )

  # Combine the plots into a grid layout
  final_plot <- ligand_plot + base_plot + receptor_plot + plot_layout(ncol = 3)

  # Save or return
  output_plot_path <- file.path(dataset_path, "figures", paste0(output_filename, ".png"))

  if (return_plot_object) {
    return(final_plot)
  } else {
    ggsave(output_plot_path, plot = final_plot, width = 16, height = 8, dpi = 300)
  }
}

plotDecipherPrioritizedMap_v4 <- function(dataset_path, receiver_cell_type = NULL, output_filename, log_transform, slice_n = 5, return_plot_object = FALSE) {
  library(ggplot2)
  library(patchwork)  # For multi-plot alignment
  
  decipher_path <- file.path(dataset_path, "data")

  # Read data ----
  lr_markers_by_cluster <- readRDS(file.path(decipher_path, "lr_markers_by_cluster.rds"))
  feature_statistics <- readRDS(file.path(decipher_path, "feature_statistics.rds"))
  decipher_scores_by_cluster <- readRDS(file.path(decipher_path, "decipher_scores_by_cluster.rds"))
  L_set <- readRDS(file.path(decipher_path, "L_set.rds"))

  ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)
  decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)

  # Normalize feature statistics ----
  normalized_feature_statistics <- feature_statistics %>%
    mutate(normalized.counts = sum.counts / n.cell) %>%
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

  if (is.null(receiver_cell_type)) {
    decipher_top_interactions_by_receiver <- decipher_scores_by_cluster_bound_clean %>%
      group_by(decipher_score_sign) %>%
      arrange(desc(abs(prioritization_score))) %>%
      select(interaction, decipher_score_sign) %>%
      distinct() %>%
      dplyr::slice_head(n = slice_n) %>%
      ungroup() %>%
      left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction"))
  } else {
    decipher_top_interactions_by_receiver <- decipher_scores_by_cluster_bound_clean %>%
      filter(receiver == receiver_cell_type) %>%
      group_by(decipher_score_sign) %>%
      arrange(desc(abs(prioritization_score))) %>%
      select(interaction, decipher_score_sign) %>%
      distinct() %>%
      dplyr::slice_head(n = slice_n) %>%
      ungroup() %>%
      left_join(decipher_scores_by_cluster_bound_clean, by = c("interaction"))
  }

  top_interactions <- decipher_top_interactions_by_receiver %>%
    select(interaction) %>%
    distinct() %>%
    unlist(use.names = FALSE)

  # Define color scale
  color_scale <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

  # Ligand (Sender) Expression Plot
  ligand_plot <- ggplot(decipher_scores_by_cluster_bound_enriched, aes(x = receiver_cluster, y = interaction)) +
    geom_point(aes(size = abs(ligand.diff.expr), fill = ligand.diff.expr), shape = 21, color = "black", stroke = 0.3) +
    color_scale +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sender Cluster", y = "Interaction", title = "Ligand Expression") +
    guides(fill = guide_colorbar(order = 1), size = guide_legend(order = 2))

  # Decipher Score Plot (Center)
  decipher_plot <- ggplot(decipher_scores_by_cluster_bound_enriched, aes(x = receiver_cluster, y = interaction)) +
    geom_point(aes(size = abs(decipher_score), fill = decipher_score), shape = 21, color = "black", stroke = 0.3) +
    color_scale +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_blank(), axis.text.y = element_blank()) +
    labs(x = "", y = "", title = "Decipher Score") +
    guides(fill = "none", size = "none")  # Hide redundant legends

  # Receptor (Receiver) Expression Plot
  receptor_plot <- ggplot(decipher_scores_by_cluster_bound_enriched, aes(x = receiver_cluster, y = interaction)) +
    geom_point(aes(size = abs(receptor.diff.expr), fill = receptor.diff.expr), shape = 21, color = "black", stroke = 0.3) +
    color_scale +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_blank(), axis.text.y = element_blank()) +
    labs(x = "Receiver Cluster", y = "", title = "Receptor Expression") +
    guides(fill = "none", size = "none")  # Hide redundant legends

  # Combine the plots into a single layout
  final_plot <- ligand_plot + decipher_plot + receptor_plot + 
    plot_layout(ncol = 3, guides = "collect") & 
    theme(legend.position = "bottom", legend.box = "vertical")  # Stack color and size legends

  # Save or return
  output_plot_path <- file.path(dataset_path, "figures", paste0(output_filename, ".png"))

  if (return_plot_object) {
    return(final_plot)
  } else {
    ggsave(output_plot_path, plot = final_plot, width = 16, height = 8, dpi = 300)
  }
}


#' Generate heatmap-ready data for a single cell type
#'
#' Constructs a long-format data frame of \`deltaPagoda\` values for a given
#' cell type across multiple conditions.  It gathers all unique regulons,
#' filters out any matching "sample", and retrieves the \`deltaPagoda\` via
#' \`get_deltaPagoda()\`.  Missing data emits warnings and yields \`NA\` values.
#'
#' @param selected_ct Character. Name of the cell type (identity) to process.
#' @param regulon_deltas_list A named list of regulon lists, one per condition.
#'   Each element should itself be a named list of data frames (by cell type),
#'   as returned by \`load_regulon_data()\`.
#' @param conditions Named list or vector whose names correspond to
#'   the conditions in \`regulon_deltas_list\`.  Used to iterate in order.
#'
#' @return A data frame with columns:
#'   - \`TF\`: regulon names  
#'   - \`Comparison\`: condition name  
#'   - \`DeltaPagoda\`: numeric \`deltaPagoda\` values (or \`NA\` if missing)
#'
#'
#' @export
generate_heatmap_data_for_celltype <- function(selected_ct, regulon_deltas_list, conditions) {
  heatmap_data <- data.frame()

  # Gather all unique regulons for this cell type across all conditions
  all_regulons <- character(0)
  for (cond_name in names(conditions)) {
    if (!is.null(regulon_deltas_list[[cond_name]]) && !is.null(regulon_deltas_list[[cond_name]][[selected_ct]])) {
      current_regulons <- regulon_deltas_list[[cond_name]][[selected_ct]]$name
      all_regulons <- union(all_regulons, current_regulons)
    } else {
        warning("No data found for cell type '", selected_ct, "' in condition '", cond_name, "'")
    }
  }
  # Remove potential "sample*" regulons if they exist
  all_regulons <- unique(all_regulons[!grepl("sample", all_regulons)])

  if(length(all_regulons) == 0) {
      warning("No valid regulons found for cell type: ", selected_ct)
      return(data.frame()) # Return empty frame
  }

  # Populate the data frame with deltaPagoda values for all regulons and conditions
  for (cond_name in names(conditions)) {
    regulon_deltas <- regulon_deltas_list[[cond_name]]
    for (regulon in all_regulons) {
      mean_delta <- get_deltaPagoda(selected_ct, regulon_deltas, regulon)
      heatmap_data <- rbind(heatmap_data, data.frame(
        TF = regulon,
        Comparison = cond_name,
        DeltaPagoda = mean_delta,
        stringsAsFactors = FALSE # Important for factor handling later
      ))
    }
  }
  return(heatmap_data)
}

#' Generate a heatmap plot of deltaPagoda values
#'
#' Creates a ggplot2 heatmap for a specified set of regulons (`tfs`) across
#' multiple conditions. The fill scale is centered at zero and bounded by
#' `±absolute_max`. Returns `NULL` if no regulons are provided.
#'
#' @param data A data frame containing at least the columns
#'   \`"TF"\`, \`"Comparison"\`, and \`"DeltaPagoda"\`.
#' @param tfs Character vector. The subset of regulon names to include in the plot.
#' @param condition_names Character vector. The factor levels and order for the
#'   \`Comparison\` axis.
#' @param absolute_max Numeric. The maximum absolute value for the fill scale;
#'   the gradient limits will be \`c(-absolute_max, absolute_max)\`.
#' @param condition_label Character. A label (e.g. "moderate" or "severe") used
#'   in the plot title to indicate by which condition the top regulons were sorted.
#' @param top_n Integer. The number of top regulons (by absolute value) used, displayed in the title.
#'
#' @return A `ggplot` object representing the heatmap, or `NULL` if `tfs` is empty.
#'
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal theme element_text ggtitle
#' @export
generate_heatmap_plot <- function(data, tfs, condition_names, absolute_max, condition_label, top_n) {
  if (length(tfs) == 0) return(NULL)

  plot_data <- data %>%
    filter(TF %in% tfs) %>%
    mutate(
      TF = factor(TF, levels = tfs),
      Comparison = factor(Comparison, levels = condition_names)
    )

  ggplot(plot_data, aes(x = Comparison, y = TF, fill = DeltaPagoda)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      na.value = "grey80", name = "TF Activity\nDelta",
      limits = c(-absolute_max, absolute_max)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = rel(1.1)),
      axis.text.y = element_text(size = rel(0.9)),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5)
    ) +
    ggtitle(paste("Top", top_n, "Regulons (Sorted by", condition_label, ")"))
}


#' Generate sorted heatmap plots for multiple cell types
#'
#' For each selected cell type, prepares heatmap data, identifies the top
#' regulons for "moderate" and "severe" conditions, and creates a pair of
#' heatmap plots (one per condition). Returns a nested list of ggplot objects.
#'
#' @param selected_receiver_cells Character vector. Names of the cell types
#'   (identities) for which to generate plots.
#' @param regulon_deltas_list Named list of regulon data per condition.
#'   Each element should itself be a named list of data frames (by cell type),
#'   as returned by `load_regulon_data()`.
#' @param conditions Named vector or list whose names correspond to the
#'   conditions (e.g., `c(moderate = NULL, severe = NULL)`); only the names
#'   are used for ordering.
#' @param top_n Integer. Number of top regulons (by absolute deltaPagoda)
#'   to include in each heatmap.
#' @param absolute_max Numeric. Maximum absolute deltaPagoda for the color
#'   scale limits in the heatmap.
#'
#' @return A nested list of ggplot objects. The top level is indexed by
#'   cell type; each sub‑list has elements `moderate_sorted` and
#'   `severe_sorted` containing the corresponding heatmap (or `NULL` if
#'   no data).
#'
#'
#' @export
generate_sorted_plots <- function(selected_receiver_cells, regulon_deltas_list, conditions, top_n, absolute_max) {
  plots <- list()
  condition_names <- names(conditions)

  for (selected_ct in selected_receiver_cells) {
    cat("Generating plots for cell type:", selected_ct, "\n")
    heatmap_data_full <- generate_heatmap_data_for_celltype(selected_ct, regulon_deltas_list, conditions)

    if (nrow(heatmap_data_full) == 0) {
      warning("Skipping plots for ", selected_ct, " due to lack of data.")
      next
    }

    plots[[selected_ct]] <- list()

    for (cond in c("mild", "severe")) {
      top_tfs <- get_top_tfs(heatmap_data_full, cond, top_n)

      plot <- generate_heatmap_plot(
        data = heatmap_data_full,
        tfs = top_tfs,
        condition_names = condition_names,
        absolute_max = absolute_max,
        condition_label = cond,
        top_n = top_n
      )

      if (is.null(plot)) {
        warning("No regulons passed filtering for '", cond, "_sorted' plot in ", selected_ct)
      }

      plots[[selected_ct]][[paste0(cond, "_sorted")]] <- plot
    }
  }

  return(plots)
}

#' Combine moderate and severe heatmaps for each cell type
#'
#' Takes a nested list of heatmap plots (as returned by `generate_sorted_plots()`)
#' and for each specified cell type, stacks a title above the side-by-side
#' “moderate” and “severe” heatmaps.  If either plot is missing, a placeholder
#' “Data Not Available” plot is used.
#'
#' @param plots A nested list indexed by cell type, each containing elements
#'   `"moderate_sorted"` and `"severe_sorted"` as ggplot objects or `NULL`.
#' @param selected_receiver_cells Character vector of cell type names to process.
#'
#' @return A named list of combined ggplot objects, one per cell type.  Each
#'   element is a patchwork object with the cell type title above the two
#'   heatmaps.
#'
#' @examples
#' \dontrun{
#' # Suppose `plots` comes from generate_sorted_plots for B_cell and T_cell
#' combined <- create_combined_plots_per_celltype(
#'   plots, 
#'   selected_receiver_cells = c("B_cell", "T_cell")
#' )
#' # Display the combined plot for B_cell
#' combined$B_cell
#' }
#'
#' @importFrom ggplot2 ggplot theme_void ggtitle theme element_text
#' @importFrom patchwork wrap_plots
#' @export
create_combined_plots_per_celltype <- function(plots, selected_receiver_cells) {
  combined_plots <- list()
  for (selected_ct in selected_receiver_cells) {
    # Check if both plots exist for this cell type
    plot_moderate <- plots[[selected_ct]][["mild_sorted"]]
    plotsevere <- plots[[selected_ct]][["severe_sorted"]]

    # Create placeholder plots if one or both are missing
    placeholder_plot <- ggplot() + theme_void() + ggtitle("Data Not Available") + theme(plot.title = element_text(hjust = 0.5))
    if (is.null(plot_moderate)) plot_moderate <- placeholder_plot
    if (is.null(plotsevere)) plotsevere <- placeholder_plot

    # Create title plot
    title <- ggplot() + theme_void() + ggtitle(selected_ct) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

    # Combine the two heatmaps side-by-side
    heatmaps <- wrap_plots(plot_moderate, plotsevere, ncol = 2)

    # Stack title above heatmaps
    combined_plots[[selected_ct]] <- wrap_plots(title, heatmaps, ncol = 1, heights = c(0.1, 1)) # Adjust height ratio for title
  }
  return(combined_plots)
}


#' Save grouped combined heatmap plots to files
#'
#' Splits a list of combined patchwork plots into groups of a specified size
#' and saves each group as a PNG file under a constructed output directory.
#'
#' @param combined_plots Named list of patchwork (ggplot) objects, one per cell type.
#' @param clusters_per_group Integer. Number of cell‑type plots to include per output file.
#' @param output_dir_base Character. Base directory path under which the output folder
#'   will be created.
#' @param output_folder_name Character. Name of the subfolder (under `output_dir_base`)
#'   where image files will be saved.
#'
#' @return Invisibly returns `NULL`. Side effects: creates the output directory (if needed)
#'   and writes PNG files named `Combined_Sorted_Heatmaps_Group_<n>.png`.
#'
#'
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom ggplot2 theme element_text unit
#' @importFrom grDevices png
#' @export
save_grouped_plots <- function(combined_plots, clusters_per_group, output_dir_base, output_folder_name) {
  # Construct the full output path
  output_dir <- file.path(output_dir_base, output_folder_name) # Match structure from file 1
  if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
  }

  cluster_names <- names(combined_plots)
  if (length(cluster_names) == 0) {
      warning("No plots available to save.")
      return()
  }

  for (i in seq(1, length(cluster_names), by = clusters_per_group)) {
    # Select clusters for the current group
    selected_clusters_group <- cluster_names[i:min(i + clusters_per_group - 1, length(cluster_names))]

    # Combine the plots for the selected group (e.g., arrange in 2 columns)
    # Adjust ncol based on how many plots per page you want vertically vs horizontally
    plots_to_save_list <- combined_plots[selected_clusters_group]

    # Check if list is empty before plotting
    if(length(plots_to_save_list) > 0){
        to_plot <- wrap_plots(plots_to_save_list, ncol = 1) + # Stack cell types vertically by default
            plot_layout(guides = "collect") & # Collect legends
            theme(
                legend.position = "bottom",
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14, face = "bold"),
                legend.key.size = unit(1.2, "cm"), # Adjust key size
                # Adjust key height/width if needed, size often covers it
                # legend.key.height = unit(1.0, "cm"),
                # legend.key.width = unit(1.0, "cm")
            )

        # Define filename
        group_num <- (i - 1) %/% clusters_per_group + 1
        file_name <- paste0("Combined_Sorted_Heatmaps_Group_", group_num, ".png")
        file_path <- file.path(output_dir, file_name)

        # Calculate dynamic height (adjust base height and multiplier as needed)
        plot_height <- 7 * length(selected_clusters_group) # Base height * number of plots stacked
        plot_width <- 12 # Fixed width (adjust if needed)

        cat("Saving group", group_num, "to:", file_path, "\n")
        # Save the plot
        ggsave(
          filename = file_path,
          plot = to_plot,
          width = plot_width,
          height = plot_height,
          dpi = 300,
          limitsize = FALSE # Important for potentially large plots
        )
    } else {
        cat("Skipping save for group", (i - 1) %/% clusters_per_group + 1, "as no plots were generated for the selected clusters.\n")
    }
  }
}


#' Create a custom color scale for flagged method points
#'
#' Generates a named vector of colors for a set of methods, assigning each method
#' two colors: one “darker” version for flagged points and the original base color
#' for unflagged points. The names follow the pattern `"MethodName.TRUE"` and
#' `"MethodName.FALSE"`.
#'
#' @param methods Character vector of method names for which to create colors.
#'   These should correspond to names in the global `method_colors` object.
#'
#' @return A named character vector of colors.  For each method `m`, the vector
#'   contains two entries:
#'   - `m.TRUE`: the darker color (for flagged = TRUE)  
#'   - `m.FALSE`: the base color  (for flagged = FALSE)
#'
#'
#' @importFrom grDevices adjustcolor
#' @export
create_flag_color_scale <- function(methods) {
  # Get a base color for each method using a hue palette.
  #base_cols <- scales::hue_pal()(length(methods))
  base_cols <- method_colors
  #names(base_cols) <- methods
  names(base_cols) <- names(method_colors)
  
  # Helper to darken a given color.
  darker <- function(col) {
    adjustcolor(col, red.f = 0.6, green.f = 0.6, blue.f = 0.6)
  }
  
  # For each method, return a vector with:
  #   - A darker color (flagged)
  #   - The base color (not flagged)
  color_values <- unlist(lapply(methods, function(m) {
    c(darker(base_cols[m]), base_cols[m])
  }))
  
  # Name the colors to match the values produced by interaction(method, flagged).
  # That is, names like "MethodName.TRUE" and "MethodName.FALSE".
  names(color_values) <- unlist(lapply(methods, function(m) {
    c(paste0(m, ".TRUE"), paste0(m, ".FALSE"))
  }))
  
  return(color_values)
}

#' Plot a Regulon Network Graph
#'
#' Builds and renders an undirected igraph network of regulon→target edges,
#' coloring and sizing vertices based on their log₂FC and whether they are
#' core regulons.
#'
#' @param data    A \code{list} with components:
#'                \describe{
#'                  \item{\code{edges}}{data.frame with \code{from}, \code{to} columns.}
#'                  \item{\code{combined_data}}{data.frame with \code{tg_gene} and \code{avg_log2FC}.}
#'                }
#' @param regulons Character vector of core regulon gene names to highlight.
#'
#' @return Invisibly returns the \code{igraph} object (after plotting).
#' @export
plot_tf_tg_network <- function(data, regulons) {
  g <- graph_from_data_frame(data$edges, directed = FALSE)
  log2fc_values <- setNames(data$combined_data$avg_log2FC, data$combined_data$tg_gene)
  log2fc_colors <- generate_log2fc_colors(log2fc_values)

  vertex_attrs <- set_vertex_attributes(g, log2fc_colors, regulons)
}




#' Generate and Save a Volcano Plot of Differential TF Activity
#'
#' @description
#' This function performs a differential analysis of transcription factor (TF)
#' activity between "case" and "control" groups for a specific cell cluster.
#' It generates a volcano plot to visualize the results, saves the plot as a PNG,
#' and saves the underlying data as a CSV file.
#'
#' @details
#' The function reads pre-calculated regulon scores and delta values, performs a
#' t-test for each TF regulon, and then uses ggplot2 to create the volcano plot.
#' It assumes the presence of a helper function `do_t_test_by_feature_by_grouping_factor()`
#' in the environment.
#'
#' @param output_data_filepath A string. The file path to the directory
#'   containing the input RDS files (`decipher_seurat_lr.rds`,
#'   `regulon_deltas_by_cluster.rds`, etc.).
#' @param selected_cluster A string. The name of the cell cluster to analyze
#'   (e.g., "CD14+ Monocytes"). This name must exist as a key in the lists
#'   within the loaded RDS files.
#' @param figures_folder A string. The path to the directory where the output
#'   plot and CSV file will be saved.
#' @param p_threshold A numeric value. The p-value cutoff for determining
#'   statistical significance. Defaults to `0.01`.
#' @param fc_threshold A numeric value. The absolute delta pagoda score
#'   (fold change) cutoff for determining biological significance.
#'   Defaults to `2.883`.
#' @param output_filename A string. The base filename for the saved volcano
#'   plot PNG image and csv.
#'
#' @return A ggplot object representing the generated volcano plot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming your data and figures folders are set up
#'   volcano_plot_object <- differentialTFActivityVolcanoPlotByCluster(
#'     output_data_filepath = "results/covid_data",
#'     selected_cluster = "CD14_plus_BDCA1_plus_cells",
#'     figures_folder = "results/figures",
#'     p_threshold = 0.05,
#'     fc_threshold = 2.0
#'   )
#'   # You can then print or further modify the returned plot object
#'   print(volcano_plot_object)
#' }
differentialTFActivityVolcanoPlotByCluster <- function(
  output_data_filepath,
  selected_cluster,
  figures_folder,
  p_threshold = 0.01,
  fc_threhsold = 2.883,
  output_filename = "diff_regulong_act_volcano_plot"){
  
  #load data
  decipher_seurat_lr <- readRDS(file.path(output_data_filepath,"decipher_seurat_lr.rds"))
  regulon_deltas_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
  regulons_scores_by_clusters <- readRDS(file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))

  #data preparation
  selected_regulons_scores_by_cluster <- regulons_scores_by_clusters[[selected_cluster]]
  condition_match <- match(colnames(selected_regulons_scores_by_cluster),names(decipher_seurat_lr$condition))
  group_vector <- decipher_seurat_lr$condition[condition_match]

  regulon_deltas_selected_cluster <- regulon_deltas_by_cluster[[selected_cluster]] %>%
  filter(class == "real")

  group_vector[group_vector == "control"] <- 0
  group_vector[group_vector == "case"] <- 1
  group_factor <- factor(c(group_vector), levels = c(0, 1), labels = c("control", "case"))

  diff_regulon_scores_p_values <- do_t_test_by_feature_by_grouping_factor(selected_regulons_scores_by_cluster,group_factor)

  regulon_deltas_selected_cluster$p_value <- diff_regulon_scores_p_values[regulon_deltas_selected_cluster$name]
  regulon_deltas_selected_cluster$log_10 <-  -1*log(regulon_deltas_selected_cluster$p_value,base=10)

  ##visualization parameters ----
  abs_max_regulon_delta <- max(abs(regulon_deltas_selected_cluster$deltaPagoda))

  ##visualization ----
  p <- ggplot(regulon_deltas_selected_cluster, aes(x = deltaPagoda, y = log_10)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", colour = "grey") +
    geom_vline(xintercept = fc_threhsold, linetype = "dashed", colour = "grey") +
    geom_vline(xintercept = -fc_threhsold, linetype = "dashed", colour = "grey") +
    geom_point(shape = 16, alpha = 0.5) +
    scale_colour_manual(values = c("0" = "#808080", "1" = "#ff8080", "2" = "#8080ff", "3" = "#ff80ff")) +
    scale_x_continuous(limits = c(-1*6,abs_max_regulon_delta), breaks = seq(-6,6,2)) +
    scale_y_continuous(limits = c(0,300)) +
    xlab("delta TF activity") +
    ylab("-log10(P)") +
    theme_classic(9) +
    theme(legend.position = "none",
          axis.title.x= element_text(size = 20),
          axis.title.y= element_text(size = 20))

  p <- p + geom_text(data = subset(regulon_deltas_selected_cluster, log_10 > -log(p_threshold,base=10) & abs(deltaPagoda) > fc_threhsold), aes(label = name),
                    vjust = "inward", hjust = "inward", check_overlap = TRUE,size=4)

  output_plot_filename <- paste(output_filename,".png")
  ggsave(file.path(figures_folder,output_plot_filename), plot = p, width = 4, height = 7, dpi = 300)

  output_data_filename <- paste(output_filename,".csv")
  write.csv(
    regulon_deltas_selected_cluster,
    file = file.path(figures_folder, output_data_filename),
    row.names = TRUE
  )

  return(p)
}

#' Plot Ligand-Receptor and Transcription Factor Interaction Heatmap
#'
#' @description
#' This function analyzes decipher scores to identify and visualize the relationships
#' between top ligand-receptor (LR) interactions and significantly altered
#' transcription factor (TF) regulons for a specific cell cluster.
#'
#' It generates a heatmap of the LR-TF decipher scores and saves it as a PNG file,
#' along with a CSV file of the underlying data matrix.
#'
#' @param output_data_filepath A character string specifying the path to the directory
#'   containing the required input RDS files (e.g., "regulon_deltas_by_cluster.rds").
#' @param selected_cluster A character string with the name of the cell cluster to analyze
#'   (e.g., "Tumour_Cells").
#' @param figures_folder A character string specifying the path to the folder where
#'   the output heatmap and CSV file will be saved.
#' @param p_threshold A numeric value for the p-value threshold used to filter for
#'   significant regulons. Defaults to `0.01`.
#' @param output_name A character string for the base name of the output PNG and CSV
#'   files (without file extension). Defaults to `"lr_tf_heatmap"`.
#' @param n_interactions A numeric value for the number of top LR interactions (ranked by
#'   absolute decipher score) to include in the heatmap. Defaults to `10`.
#' @param min_abs_decipher_score A numeric value for the minimum absolute decipher score
#'   for filtering interactions. Defaults to `0.4`.
#' @param min_delta_pagoda A numeric value for the minimum absolute `deltaPagoda` score
#'   used to filter for significant regulons. Defaults to `2`.
#'
#' @return The function does not return any R object. It saves two files to the
#'   `figures_folder`:
#'   \itemize{
#'     \item A PNG file containing the LR-TF heatmap.
#'     \item A CSV file containing the matrix of decipher scores used to generate the heatmap.
#'   }
#'
#' @importFrom dplyr filter mutate rename select arrange distinct slice_head ungroup left_join bind_rows if_else pull
#' @importFrom reshape2 acast
#' @importFrom gplots heatmap.2
#' @importFrom grDevices png dev.off
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This is a hypothetical example, as it requires specific input files.
#' plotLRTFHeatmap(
#'   output_data_filepath = "path/to/your/data_output",
#'   selected_cluster = "My_Cell_Cluster",
#'   figures_folder = "path/to/your/figures",
#'   p_threshold = 0.05,
#'   n_interactions = 15
#' )
plotLRTFHeatmap <- function(
  output_data_filepath,
  selected_cluster,
  figures_folder,
  p_threshold = 0.01,
  output_name = "lr_tf_heatmap",
  n_interactions = 10,
  min_abs_decipher_score = 0.4,
  min_delta_pagoda = 2
){
  regulon_deltas_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
  regulons_scores_by_clusters <- readRDS(file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))
  decipher_seurat_lr <- readRDS(file.path(output_data_filepath,"decipher_seurat_lr.rds"))
  decipher_scores_by_cluster <- readRDS(file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
  decipher_scores_by_regulon_and_cluster <- readRDS(file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))

  ##data wrangling ----
  regulon_scores_selected_cluster <- regulons_scores_by_clusters[[selected_cluster]]
  condition_match <- match(colnames(regulon_scores_selected_cluster),names(decipher_seurat_lr$condition))
  group_vector <- decipher_seurat_lr$condition[condition_match]

  regulon_deltas_selected_cluster <- regulon_deltas_by_cluster[[selected_cluster]] %>%
    filter(class == "real")

  group_vector[group_vector == "control"] <- 0
  group_vector[group_vector == "case"] <- 1
  group_factor <- factor(c(group_vector), levels = c(0, 1), labels = c("control", "case"))

  diff_regulon_scores_p_values <- do_t_test_by_feature_by_grouping_factor(regulon_scores_selected_cluster,group_factor)

  regulon_deltas_selected_cluster$p_value <- diff_regulon_scores_p_values[regulon_deltas_selected_cluster$name]
  regulon_deltas_selected_cluster$log_10 <-  -1*log(regulon_deltas_selected_cluster$p_value,base=10)

  regulon_signature <- regulon_deltas_selected_cluster %>%
    filter(log_10 > -log(p_threshold,base=10) & abs(deltaPagoda) > min_delta_pagoda) %>%
    pull(name)

  decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)
  decipher_scores_by_cluster_bound_filtered <- decipher_scores_by_cluster_bound %>%
    mutate(decipher_score = sign(decipher_score)*log10(abs(decipher_score)+1)) %>%
    filter(abs(decipher_score) > min_abs_decipher_score)

  decipher_scores_by_cluster_bound_clean <- decipher_scores_by_cluster_bound %>%
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    rename(receiver=receiver_cluster,sender=sender_cluster,prioritization_score=decipher_score) %>%
    select(interaction,ligand,receptor,receiver,prioritization_score) %>%
    arrange(prioritization_score)


  decipher_top_interactions_cluster_selected_cluster <- decipher_scores_by_cluster_bound_clean %>%
    #filter(receiver %in% "Tumour_Cells") %>%
    mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
    #group_by(decipher_score_sign) %>%
    arrange(desc(abs(prioritization_score))) %>%
    select(interaction) %>%
    distinct() %>%
    slice_head(n = n_interactions) %>%
    ungroup() %>%
    left_join(decipher_scores_by_cluster_bound)

  top_interactions <- decipher_top_interactions_cluster_selected_cluster %>%
    select(interaction) %>%
    distinct() %>%
    unlist(use.names=FALSE)

  to_plot <- decipher_scores_by_regulon_and_cluster[[selected_cluster]] %>%
    filter(regulon %in% regulon_signature,
          interaction %in% top_interactions)


  this_matrix <-reshape2::acast(to_plot,interaction~regulon,value.var = "decipher_score",fill = 0)

  ##visualization ----
  file_name <- paste(output_name,".png",sep="")
  png(file.path(figures_folder,file_name),width = 15, height = 9, units = "cm",res=600)
  heatmap.2(
    this_matrix,
    trace="none",
    col = "bluered",
    breaks = 100,
    cexRow = 0.7,
    cexCol=0.5,
    scale = "none",
    key.title = "LR-TF Decipher score",Rowv = FALSE,
    dendrogram = "none",
    margins = c(6,6),
    key = FALSE,
    keysize = 0.3,
    Colv=TRUE)
  
  dev.off()
  file_name <- paste(output_name,".csv",sep="")
  write.csv(this_matrix,file.path(figures_folder,file_name),row.names=TRUE)
}



#' Plot cluster-sorted TF heatmap
#'
#' Select the top `top_n` transcription factors (TFs) by absolute `deltaPagoda`
#' in a `selected_cluster`, then plot a heatmap of those TFs across *all*
#' clusters in `deltas_by_cluster`. The TF ordering is taken from the selected
#' cluster and applied to every cluster.
#'
#' Each element of `deltas_by_cluster` must be a named list element (cluster name)
#' whose value is either:
#' * a `data.frame` with columns `name` and `deltaPagoda`, or
#' * a named numeric vector (names = TF names, values = deltas).
#'
#' NA values are allowed and shown as `grey80` in the plot. If `global_max` is
#' provided it fixes the symmetric color scale to `[-global_max, +global_max]`.
#'
#' @param selected_cluster Character scalar. Name of the cluster to use for
#'   selecting and ordering the top TFs.
#' @param top_n Integer. Number of top TFs (by absolute delta in
#'   `selected_cluster`) to select. Default 20.
#' @param deltas_by_cluster Named list. Each element is a cluster (see details).
#' @param global_max Numeric scalar or `NULL`. If numeric, forces symmetric
#'   color limits to `c(-global_max, global_max)`. If `NULL` (default) limits
#'   are taken from the data.
#' @param cluster_order Optional character vector specifying desired cluster row
#'   order. Names not present in `deltas_by_cluster` are warned and ignored; any
#'   remaining clusters are appended after the provided ordering.
#'
#' @return A `ggplot2` object (tile heatmap) with TFs on the x-axis and Clusters
#'   on the y-axis.
#'
#' @examples
#' \dontrun{
#' my_deltas_list <- list(
#'   CM_CD8  = data.frame(name = c("TF1","TF2","TF3"), deltaPagoda = c(0.5, -1.2, 0.2)),
#'   EM_CD8  = c(TF1 = 0.3, TF2 = -0.8, TF3 = 0.1),
#'   GZMK_CD8 = data.frame(name = c("TF1","TF2"), deltaPagoda = c(0.2, NA))
#' )
#'
#' p <- plotClusterSortedTFHeatmap(
#'   selected_cluster = "CM_CD8",
#'   top_n = 10,
#'   deltas_by_cluster = my_deltas_list,
#'   global_max = 1.5
#' )
#' print(p)
#' }
#'
#' @seealso \code{\link[ggplot2]{geom_tile}}, \code{\link[tibble]{enframe}}
#' @keywords visualization heatmap TF
#' @export
#' @importFrom tibble enframe
#' @importFrom dplyr bind_rows filter arrange desc slice_head pull mutate
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal theme element_text
plotClusterSortedTFHeatmap <- function(
  selected_cluster,
  top_n = 20,
  deltas_by_cluster,
  global_max = NULL,
  cluster_order = NULL   # optional character vector giving desired cluster row order
) {
  # helper: coerce element to named numeric vector
  vec_from <- function(x) {
    if (is.data.frame(x)) {
      if (all(c("deltaPagoda","name") %in% names(x))) {
        v <- x$deltaPagoda; names(v) <- x$name; return(v)
      } else stop("Data.frame must have columns 'deltaPagoda' and 'name'.")
    } else if (is.numeric(x) && !is.null(names(x))) {
      return(x)
    } else stop("Each element must be a named numeric vector or a data.frame with 'deltaPagoda'+'name'.")
  }

  if (!is.list(deltas_by_cluster) || length(deltas_by_cluster) == 0) stop("deltas_by_cluster must be a non-empty named list of clusters.")
  if (is.null(names(deltas_by_cluster)) || any(names(deltas_by_cluster) == "")) stop("deltas_by_cluster must be a named list (cluster names).")
  if (!selected_cluster %in% names(deltas_by_cluster)) stop("selected_cluster not found in deltas_by_cluster.")

  # build long table: TF, Delta, Cluster
  cluster_names <- names(deltas_by_cluster)
  rows <- lapply(cluster_names, function(cl) {
    v <- tryCatch(vec_from(deltas_by_cluster[[cl]]), error = function(e) {
      warning("Skipping cluster ", cl, ": ", e$message); return(NULL)
    })
    if (is.null(v)) return(tibble::tibble(TF = character(), Delta = numeric(), Cluster = character()))
    tibble::enframe(v, name = "TF", value = "Delta") %>% dplyr::mutate(Cluster = cl)
  })
  df_long <- dplyr::bind_rows(rows)

  # pick top TFs by |Delta| in selected cluster
  top_tfs <- df_long %>%
    dplyr::filter(Cluster == selected_cluster, !is.na(Delta)) %>%
    dplyr::arrange(dplyr::desc(abs(Delta))) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(Delta) %>%    # order so negatives are left, positives right
    dplyr::pull(TF)

  if (length(top_tfs) == 0) stop("No TFs found for selected cluster (or all NAs).")

  # subset & set factor levels: TFs by top_tfs; clusters by cluster_order or as present
  if (!is.null(cluster_order)) {
    missing <- setdiff(cluster_order, cluster_names)
    if (length(missing)) warning("cluster_order contains names not in deltas_by_cluster: ", paste(missing, collapse=", "))
    cluster_levels <- intersect(cluster_order, cluster_names)
    cluster_levels <- c(cluster_levels, setdiff(cluster_names, cluster_levels)) # keep rest after
  } else {
    cluster_levels <- cluster_names
  }

  df_top <- df_long %>%
    dplyr::filter(TF %in% top_tfs) %>%
    dplyr::mutate(
      TF = factor(TF, levels = top_tfs),
      Cluster = factor(Cluster, levels = cluster_levels)
    )

  # color limits
  if (!is.null(global_max) && is.numeric(global_max) && length(global_max) == 1) {
    limits <- c(-abs(global_max), abs(global_max))
  } else {
    lim <- max(abs(df_top$Delta), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 1
    limits <- c(-lim, lim)
  }

  # plot
  p <- ggplot2::ggplot(df_top, ggplot2::aes(x = TF, y = Cluster, fill = Delta)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                                  limits = limits, na.value = "grey80") +
    ggplot2::labs(title = paste0(selected_cluster, " top ", length(top_tfs), " TFs (sorted by ", selected_cluster, ")"),
                  x = "TF", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = element_blank(),
                   legend.position = "bottom")

  return(p)
}

#' Plot PubMed-prioritized TF→TG heatmaps per cluster
#'
#' For each cluster in `selected_clusters`, this function:
#' * reads precomputed objects from `output_data_filepath` (`capped_regulons_all_clusters.rds`,
#'   `regulon_deltas_by_cluster.rds`, `significant_regulon_markers_by_cluster.rds`),
#' * selects the top regulons (class == "real", ordered by `deltaPagoda`),
#' * extracts regulon data via `extract_regulon_data()` and (optionally) calls
#'   `plot_tf_tg_network()` for a network visualization,
#' * builds a binary TF → TG matrix, filters target genes with >1 incoming TF,
#' * queries PubMed counts via `get_n_pubmed_articles_per_gene()` and selects
#'   the top `n_pubmed` genes,
#' * marks those top genes in the matrix (1 → 1.2), clusters rows/columns,
#'   plots a heatmap and saves a PNG named `<cluster>_tgs_heatmap_pubmed.png`
#'   into `figures_folder`.
#'
#' **Important:** this function calls external helper functions that must be
#' available in the environment:
#' `extract_regulon_data`, `plot_tf_tg_network`, and
#' `get_n_pubmed_articles_per_gene`. It also expects a variable
#' `figures_folder` to be defined (or defined in the calling environment).
#'
#' @param output_data_filepath Character scalar. Path to a directory containing
#'   the RDS files:
#'   - `capped_regulons_all_clusters.rds`
#'   - `regulon_deltas_by_cluster.rds`
#'   - `significant_regulon_markers_by_cluster.rds`
#' @param selected_clusters Character vector. Names of clusters to process (these
#'   must be keys in the `regulon_deltas_by_cluster` object).
#' @param n_pubmed Integer (default 40). Number of top genes by PubMed mention
#'   frequency to keep for the heatmap.
#'
#' @return Invisibly returns \code{NULL}. Primary purpose is side-effects:
#'   saving PNG heatmaps to the \code{figures_folder} path and printing messages.
#'
#' @details
#' - NA values in PubMed counts or missing edges may lead to clusters being
#'   skipped with a warning.
#' - The function uses \code{reshape2::melt()} in the body; if you prefer the
#'   tidyverse variant, replace \code{reshape2::melt()} by
#'   \code{tidyr::pivot_longer()}.
#' - The plot color scale uses a three-color gradient (white → black → red)
#'   and the values are slightly bumped (1 -> 1.2) to emphasize top genes.
#'
#' @examples
#' \dontrun{
#' # define the folder where figures should be written (must exist)
#' figures_folder <- "figures"
#'
#' # run the function (requires helper functions to be defined)
#' plot_pubmed_tg_heatmaps(
#'   output_data_filepath = "results/2025-09-25",
#'   selected_clusters = c("CM_CD8", "EM_CD8"),
#'   n_pubmed = 40
#' )
#' }
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[ggplot2]{geom_tile}},
#'   \code{\link[tidyr]{pivot_wider}}
#' @keywords visualization heatmap pubmed transcription-factors
#' @export
plot_pubmed_tg_heatmaps <- function(
  output_data_filepath,
  selected_clusters,
  n_pubmed = 40
){
  capped_regulons_all_clusters <- readRDS(file.path(output_data_filepath,"capped_regulons_all_clusters.rds"))
  regulon_deltas_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
  significant_regulon_markers_by_cluster <- readRDS(file.path(output_data_filepath,"significant_regulon_markers_by_cluster.rds"))
  for(selected_ct in selected_clusters){
    # Run the refactored code
    regulons <- regulon_deltas_by_cluster[[selected_ct]] %>%
      filter(class == "real") %>%
      arrange(desc(deltaPagoda)) %>%
      dplyr::slice_head(n=10) %>%
      pull(name)
    data <- extract_regulon_data(regulons, capped_regulons_all_clusters, significant_regulon_markers_by_cluster, selected_ct)
    #plot_tf_tg_network(data, regulons)

    # Create a binary matrix
    binary_matrix <- data$edges %>%
      mutate(value = 1) %>%  # Assign 1 for each edge
      tidyr::pivot_wider(names_from = to, values_from = value, values_fill = 0) %>%
      tibble::column_to_rownames("from")

    # Filter columns where more than one 'from' talks to a 'to'
    binary_matrix <- binary_matrix[, colSums(binary_matrix) > 1]

    #see frequency of genes mentioned in pubmed
    # Define genes
    genes <- colnames(binary_matrix)  # Replace with your genes of interest

    # # Apply the function to each gene and store the results in a data frame
    PubMed_Count = sapply(genes, get_n_pubmed_articles_per_gene)

    gene_counts <- data.frame(
      Gene = genes,
      PubMed_Count = PubMed_Count
    )

    top_genes <- gene_counts %>%
      filter(Gene %in% genes) %>%
      arrange(desc(PubMed_Count)) %>%
      top_n(n = n_pubmed, wt = PubMed_Count) %>%
      pull(Gene)

    binary_matrix <- binary_matrix[,top_genes]

    # Convert binary_matrix to matrix format
    binary_matrix_mat <- as.matrix(binary_matrix)
    # Define the color palette with custom logic for red and white
    for(i in 1:nrow(binary_matrix_mat)){
      for(j in 1:ncol(binary_matrix_mat)){
        if(binary_matrix_mat[i,j] == 1 & colnames(binary_matrix_mat)[j] %in% top_genes)
          binary_matrix_mat[i,j] <- 1.2
      }
    }
    my_palette <- c("white", "black", "red")

    # Perform hierarchical clustering on rows and columns
    row_dendrogram <- hclust(dist(binary_matrix_mat))   # Hierarchical clustering on rows
    col_dendrogram <- hclust(dist(t(binary_matrix_mat))) # Hierarchical clustering on columns

    # Order the matrix based on clustering
    binary_matrix_ordered <- binary_matrix_mat[row_dendrogram$order, col_dendrogram$order]

    # Convert the ordered matrix to a data frame for ggplot2
    binary_matrix_df <- as.data.frame(binary_matrix_ordered)
    binary_matrix_df$row <- factor(rownames(binary_matrix_ordered), levels = rownames(binary_matrix_ordered)) # Set ordered row levels

    # Melt the data for ggplot2
    binary_matrix_melted <- reshape2::melt(binary_matrix_df, id.vars = "row")

    # Create the heatmap plot
    p <- ggplot(binary_matrix_melted, aes(x = variable, y = row, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colors = my_palette,
                          breaks = c(-0.1, 0.9, 1.1, 1.2),
                          limits = c(-0.1, 1.2),
                          guide = "none") +
      labs(title = NULL) +
      theme_minimal() +
      theme(
        #plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        # Conditionally set y-axis labels
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "lines")
      ) +
      scale_y_discrete(position = "right") #do this for severe and comment out for moderate
      #scale_y_discrete(labels = function(x) ifelse(x %in% top_genes, x, ""))  # Conditional row labels
    filename <- paste0(selected_ct,"_tgs_heatmap_pubmed.png")
    ggsave(file.path(figures_folder,filename),p,width = 12,height = 5.5,units="cm")

    return(p)
  }
}

#' Plot a layered network graph and save to PNG
#'
#' Render an `igraph` object with a simple deterministic layered layout
#' (Sender Cell Type → Ligand → Receptor → TF). Nodes within each layer are
#' evenly spaced horizontally; layers are stacked vertically. The function
#' writes a PNG to `output_file`.
#'
#' @param g An `igraph` object. Must have a vertex attribute `layer` (one of
#'   `"Sender Cell Type"`, `"Ligand"`, `"Receptor"`, `"TF"`) and vertex
#'   attributes used for plotting (e.g. `color`, `size`).
#' @param output_file Character scalar. Path to the output PNG file to write.
#' @param cluster_name Character scalar. Cluster name used for the plot title.
#' @param condition_label Character scalar. Additional label used in the title
#'   (e.g. `"case"` or `"control"`).
#'
#' @details
#' - The layout is deterministic and simple: nodes are placed in one of four
#'   horizontal bands according to `V(g)$layer`. Within each band nodes are
#'   evenly spaced left→right. If a vertex's `layer` is not one of the four
#'   expected values it will be omitted from the layout (left as `NA`).
#' - The function uses base `png()` and `plot.igraph()` and writes the file at
#'   high resolution (400 DPI). Existing files with the same name will be
#'   overwritten.
#'
#' @return Invisibly returns \code{NULL}. Primary effect is the written PNG file.
#'
#' @examples
#' \dontrun{
#' g <- create_graph(all_edges, nodes)
#' plot_graph(g, "figs/CM_CD8_case_network.png", "CM_CD8", "case")
#' }
#'
#' @keywords plotting network igraph
#' @export
#' @importFrom igraph V
plot_graph <- function(g, output_file, cluster_name, condition_label) {
  #layout <- layout_with_sugiyama(g)$layout
  #layout <- layout_as_tree(g, root = V(g)[V(g)$layer == "Sender Cell Type"], circular = FALSE)
  # Custom layout: even horizontal spacing per layer
  node_layers <- igraph::V(g)$layer
  unique_layers <- c("Sender Cell Type", "Ligand", "Receptor", "TF")

  layout <- matrix(NA, nrow = length(igraph::V(g)), ncol = 2)

  for (i in seq_along(unique_layers)) {
    layer <- unique_layers[i]
    nodes_in_layer <- which(node_layers == layer)
    n <- length(nodes_in_layer)
    if (n > 0) {
      x_pad <- 4  # stretch width of layout
      #x_positions <- seq(from = 0, to = x_pad, length.out = n)
      #x_positions <- seq(from = 0, to = n * 0.3, length.out = n)
      x_positions <- seq(from = 0, to = 1, length.out = n + 2)[2:(n + 1)]
      y_position <- length(unique_layers) - i
      layout[nodes_in_layer, 1] <- x_positions
      layout[nodes_in_layer, 2] <- y_position
    }
  }

  png(output_file, width = 10, height = 6.7, units = "in", res = 400)
  plot(g,
       layout = layout,
       vertex.label = igraph::V(g)$name,
       vertex.label.cex = 0.9,       # Font size
       vertex.label.font = 2,        # Font weight (2 = bold)
       vertex.label.color = "black",
       main = paste(cluster_name, " (", condition_label, ")", sep = " "),
       edge.arrow.size = 0.5,
       vertex.frame.color = NA,
       asp = 0.5)
  dev.off()
  message("Saved network plot to: ", output_file)

}

#' Generate and save a Sender→Ligand→Receptor→TF network plot for a cluster/condition
#'
#' Orchestrates the full pipeline to build a communication network for a given
#' `cluster_name` and `condition_label` and writes CSVs and a PNG image to
#' `output_dir`. The function:
#' 1. validates inputs (`check_cluster_exists()`),  
#' 2. selects top regulons (`get_top_regulons()`),  
#' 3. selects top interactions (`get_top_interactions()`),  
#' 4. builds Sender→Ligand edges (`build_sender_ligand_edges()`),  
#' 5. builds Ligand→Receptor edges (`build_ligand_receptor_edges()`),  
#' 6. builds Receptor→TF edges (`build_receptor_tf_edges()`),  
#' 7. composes node/edge tables (`build_graph_components()`),  
#' 8. creates an `igraph` (`create_graph()`),  
#' 9. assigns edge colours/widths (`assign_edge_colors()`), and  
#' 10. plots and saves the graph (`plot_graph()`).
#'
#' Many helper functions are called and therefore must be present in the
#' environment: `check_cluster_exists`, `get_top_regulons`,
#' `get_top_interactions`, `extract_regulon_data`, `plot_tf_tg_network`,
#' `build_sender_ligand_edges`, `build_ligand_receptor_edges`,
#' `build_receptor_tf_edges`, `build_graph_components`, `create_graph`,
#' `assign_edge_colors`, and `plot_graph`. The function has side-effects: it
#' writes `<label>_<cluster>_edges.csv`, `<label>_<cluster>_nodes.csv` and a PNG
#' network image into `output_dir`.
#'
#' @param condition_label Character scalar. Raw condition key (e.g. `"SevCOVID_Azimuthl2"`, `"MilCOVID_Azimuthl2"`) used to derive a pretty label.
#' @param cluster_name Character scalar. Cluster key to analyse.
#' @param decipher_scores Named list of per-cluster decipher score tables.
#' @param decipher_scores_by_regulon_and_cluster Named list of per-cluster decipher tables broken out by regulon.
#' @param regulon_deltas_by_cluster Named list of per-cluster regulon delta tables (must contain `class`, `name`, `deltaPagoda`).
#' @param feature_statistics Data.frame/tibble with per-feature stats containing at least `sum.counts`, `n.cell`, `condition`, `cluster`, `feature`.
#' @param sender_cts Character vector of cluster names to consider as candidate senders.
#' @param output_dir Character scalar. Directory where CSVs and PNG will be written (created if missing).
#' @param top_interactions Optional precomputed interactions table (passed to `get_top_interactions()`); if `NULL` the top interactions are selected automatically.
#' @param global_deltaPagoda_max Numeric scalar (or a global variable) used to scale TF node colours.
#' @param global_receptor_tf_col_max Numeric scalar used to scale receptor→TF colour mapping.
#' @param global_sender_ligand_max Numeric scalar used to scale sender→ligand weights.
#' @param global_decipher_score_max Numeric scalar used to scale ligand→receptor edge weights.
#' @param n_top_regulons Integer (default 10). Number of regulons to select per cluster.
#'
#' @return Invisibly returns \code{NULL}. Main effects are written files and a saved PNG; the function also prints progress messages.
#'
#' @details
#' - `pretty_label` and `pretty_cluster` are created via a small `switch()` mapping for nicer filenames and titles. Extend those mappings if you have more keys.
#' - If any helper step fails for the cluster (missing data, empty edges), the function will error out (or the helper will throw); consider wrapping calls if you prefer graceful skipping.
#' - This function expects `get_n_pubmed_articles_per_gene` (and other domain-specific helpers) to be available in the environment used by the helper functions called inside.
#'
#' @examples
#' \dontrun{
#' generate_network_plot(
#'   condition_label = "SevCOVID_Azimuthl2",
#'   cluster_name = "CD14_Mono",
#'   decipher_scores = decipher_scores,
#'   decipher_scores_by_regulon_and_cluster = decipher_scores_by_regulon_and_cluster,
#'   regulon_deltas_by_cluster = regulon_deltas_by_cluster,
#'   feature_statistics = feature_statistics,
#'   sender_cts = c("CD14_Mono","CD16_Mono"),
#'   output_dir = "figures",
#'   n_top_regulons = 10
#' )
#' }
#'
#' @seealso check_cluster_exists, get_top_regulons, get_top_interactions,
#'   build_sender_ligand_edges, build_ligand_receptor_edges,
#'   build_receptor_tf_edges, build_graph_components, create_graph,
#'   assign_edge_colors, plot_graph
#' @keywords network visualization ligand receptor TF
#' @export
generate_network_plot <- function(condition_label, cluster_name, 
                                  output_data_filepath,
                                  sender_cts, output_dir,top_interactions = NULL,
                                  scaling_global_deltaPagoda_max = 1,
                                  scaling_global_receptor_tf_col_max = 1,
                                  scaling_global_sender_ligand_max = 1,
                                  scaling_global_decipher_score_max = 1,
                                  n_top_regulons = 10) {

  #0. read data and calculate some plot parameters
  decipher_scores <- readRDS(file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
  decipher_scores_by_regulon_and_cluster <- readRDS(file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
  regulon_deltas_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
  feature_statistics <- readRDS(file.path(output_data_filepath,"feature_statistics.rds"))

  global_deltaPagoda_max <- max(sapply(regulon_deltas_by_cluster[target_clusters], function(x) max(x$deltaPagoda, na.rm = TRUE)))
  global_deltaPagoda_max = global_deltaPagoda_max*scaling_global_deltaPagoda_max
  # Receptor→TF imp.perm * sign(spearman.cor)
  global_receptor_tf_col_max <- max(
    sapply(c(decipher_scores_by_regulon_and_cluster[target_clusters]), function(cluster_df) {
      if (!is.null(cluster_df)) {
        df <- cluster_df %>% mutate(col = imp.perm * sign(spearman.cor))
        max(abs(df$col), na.rm = TRUE)
      } else {
        0
      }
    })
  )
  global_receptor_tf_col_max=global_receptor_tf_col_max*scaling_global_receptor_tf_col_max
  # Sender→Ligand frac.normalized.counts

  global_sender_ligand_max <- max(
    feature_statistics %>% filter(cluster %in% target_clusters) %>% pull(sum.counts) / feature_statistics %>% filter(cluster %in% target_clusters) %>% pull(n.cell),
    na.rm = TRUE
  )
  global_sender_ligand_max <- 1
  global_sender_ligand_max=global_sender_ligand_max*scaling_global_sender_ligand_max
  # Ligand→Receptor decipher_score
  global_decipher_score_max <- max(
    sapply(decipher_scores[target_clusters], function(x) max(abs(x$decipher_score), na.rm = TRUE)))
  global_decipher_score_max=global_decipher_score_max*scaling_global_decipher_score_max

  # 1. Check inputs
  check_cluster_exists(cluster_name, decipher_scores, decipher_scores_by_regulon_and_cluster, regulon_deltas_by_cluster)
  
  # 2. Get top regulons
  top_regulons_out <- get_top_regulons(cluster_name, regulon_deltas_by_cluster, n = n_top_regulons)
  top_regulons_df <- top_regulons_out$top_regulons_df
  top_regulons <- top_regulons_out$top_regulons
  
  # 3. Get top interactions
  top_interactions <- get_top_interactions(cluster_name, decipher_scores, decipher_scores_by_regulon_and_cluster, top_regulons,top_interactions = top_interactions)
  
  # 4. Build Sender → Ligand edges
  sl_out <- build_sender_ligand_edges(top_interactions, feature_statistics, sender_cts,global_sender_ligand_max)
  sender_ligand <- sl_out$sender_ligand
  edges_sender_ligand <- sl_out$edges_sender_ligand
  
  # 5. Build Ligand → Receptor edges
  this_decipher_scores <- decipher_scores[[cluster_name]]
  ligand_receptor_edges <- build_ligand_receptor_edges(this_decipher_scores, sender_ligand, top_interactions,global_decipher_score_max)
  
  # 6. Build Receptor → TF edges
  receptor_tf_edges <- build_receptor_tf_edges(top_interactions,global_receptor_tf_col_max)
  
  # 7. Combine edges and build nodes
  graph_components <- build_graph_components(edges_sender_ligand, ligand_receptor_edges, receptor_tf_edges, top_interactions, sender_ligand, top_regulons_df,global_deltaPagoda_max)
  all_edges <- graph_components$all_edges
  nodes <- graph_components$nodes

  write.csv(
    all_edges,
    file.path(output_dir,
              paste0(condition_label, "_", cluster_name, "_edges.csv")),
    row.names = TRUE
  )
  write.csv(
    nodes,
    file.path(output_dir,
              paste0(condition_label, "_", cluster_name, "_nodes.csv")),
    row.names = TRUE
  )
  
  # 8. Create graph
  g <- create_graph(all_edges, nodes)
  
  # 9. Assign edge colors and widths
  g <- assign_edge_colors(g, all_edges,global_receptor_tf_col_max)
  
  # 10. Plot and save the network
  output_file <- file.path(output_dir, paste0(condition_label, "_", cluster_name, "_network_map.png"))
  plot_graph(g, output_file, cluster_name, condition_label)
}