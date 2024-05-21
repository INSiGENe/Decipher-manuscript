# working with cross species functions
#' Convert Human Ligand-Receptor Set to Mouse
#'
#' This function converts a human ligand-receptor set to its mouse equivalent.
#' It utilizes ortholog mapping to find corresponding mouse genes for human ligands and receptors.
#' The function requires an input data frame with 'ligand' and 'receptor' columns.
#'
#' @param human_l_set A data frame representing the human ligand-receptor set.
#'        This data frame must contain the columns 'ligand' and 'receptor'.
#'
#' @return A data frame representing the mouse ligand-receptor set.
#'         The output includes 'ligand' and 'receptor' columns for mouse genes,
#'         along with 'old_ligand_name' and 'old_receptor_name' columns to reference
#'         the original human genes. It also includes an 'interaction' column
#'         representing the ligand-receptor pairs.
#'
#' @importFrom dplyr left_join rename select distinct mutate
#' @importFrom babelgene orthologs
#'
#' @export
#' @examples
#' # Example usage:
#' human_l_set <- data.frame(
#'   ligand = c("LIG1", "LIG2"),
#'   receptor = c("REC1", "REC2")
#' )
#' mouse_l_set <- convertLsetToMouse(human_l_set)
#'
convertLsetToMouse <- function(human_l_set) {
  # Ensure the input data frame has the required columns
  required_columns <- c("ligand", "receptor")
  if (!all(required_columns %in% names(human_l_set))) {
    stop("The input data frame must contain columns 'ligand' and 'receptor'.")
  }

  # Obtain mouse orthologs for the unique human ligands and receptors
  unique_genes <- unique(c(human_l_set$ligand, human_l_set$receptor))
  mouse_Lset_genes <- babelgene::orthologs(genes = unique_genes, species = "mouse")

  # Create a mapping table from human symbols to mouse symbols
  mapping_table <- data.frame(
    human_symbol = mouse_Lset_genes$human_symbol,
    mouse_symbol = mouse_Lset_genes$symbol
  )

  # Convert the human ligand-receptor set to mouse
  Lset_mouse <- human_l_set %>%
    dplyr::left_join(mapping_table, by = c("ligand" = "human_symbol")) %>%
    dplyr::rename(old_ligand_name = ligand, ligand = mouse_symbol) %>%
    dplyr::left_join(mapping_table, by = c("receptor" = "human_symbol")) %>%
    dplyr::rename(old_receptor_name = receptor, receptor = mouse_symbol) %>%
    dplyr::select(ligand, receptor, old_ligand_name, old_receptor_name) %>%
    dplyr::mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    dplyr::distinct()

  return(Lset_mouse)
}


#' Convert Human Gene Symbols to Mouse Orthologs
#'
#' This function takes a vector of human gene symbols and finds their corresponding
#' mouse orthologs using the babelgene package. It returns a vector of mouse gene
#' symbols that are orthologs of the given human genes.
#'
#' @param human_symbols A vector of human gene symbols to be converted.
#'
#' @return A vector of mouse gene symbols corresponding to the human input.
#'         NA values and duplicates are removed from the output.
#'
#' @importFrom babelgene orthologs
#' @importFrom stats na.omit
#' @examples
#' # Example usage:
#' human_genes <- c("BRCA1", "TP53")
#' mouse_genes <- convertHumanSymbolsToMouse(human_genes)
#'
#' @export
convertHumanSymbolsToMouse <- function(human_symbols) {
  # Validate input
  if (!is.vector(human_symbols) || is.null(human_symbols)) {
    stop("human_symbols must be a non-null vector.")
  }

  # Retrieve mouse orthologs for the given human gene symbols
  mouse_orthologs <- babelgene::orthologs(genes = human_symbols, species = "mouse")

  # Create a mapping table from human symbols to mouse symbols
  mapping_table <- data.frame(
    human_symbol = mouse_orthologs$human_symbol,
    mouse_symbol = mouse_orthologs$symbol
  )

  # Map human symbols to mouse symbols
  # and handle cases where there are no matching mouse symbols
  mapped_mouse_symbols <- mapping_table$mouse_symbol[match(human_symbols, mapping_table$human_symbol)]

  # Remove NA values and duplicates
  mouse_symbols <- na.omit(mapped_mouse_symbols)
  mouse_symbols <- unique(mouse_symbols)

  return(mouse_symbols)
}

#' Clean Symbols in a String
#'
#' This function replaces or removes various symbols in a given string with safer alternatives.
#' It is designed to handle common symbols that might cause issues in text processing or analysis.
#'
#' @param string A character string that needs to be cleaned.
#' @return A cleaned character string with special symbols replaced or removed.
#' @importFrom stringr str_replace_all str_remove_all
#' @export
#' @examples
#' cleanSymbols("Example string with symbols like %, &, $, etc.")
#' # [1] "Example_string_with_symbols_like_percent_and_dollar_etc"
#'
cleanSymbols <- function(string) {
  # Remove or replace various symbols with safer alternatives
  string <- stringr::str_remove_all(string, "/")
  string <- stringr::str_replace_all(string, " ", "_")
  string <- stringr::str_replace_all(string, "\\+", "_plus_")
  string <- stringr::str_replace_all(string, "\\-", "_minus_")
  string <- stringr::str_replace_all(string, "\\(", "")
  string <- stringr::str_replace_all(string, "\\)", "")
  string <- stringr::str_replace_all(string, "%", "_percent_")
  string <- stringr::str_replace_all(string, "\\.", "_dot_")
  string <- stringr::str_replace_all(string, ",", "_comma_")
  string <- stringr::str_replace_all(string, ":", "_colon_")
  string <- stringr::str_replace_all(string, ";", "_semicolon_")
  string <- stringr::str_replace_all(string, "&", "_and_")
  string <- stringr::str_replace_all(string, "\\?", "_question_")
  string <- stringr::str_replace_all(string, "!", "_exclamation_")
  string <- stringr::str_replace_all(string, "\"", "_quote_")
  string <- stringr::str_replace_all(string, "'", "_apostrophe_")
  string <- stringr::str_replace_all(string, "=", "_equals_")
  string <- stringr::str_replace_all(string, "\\*", "_asterisk_")
  string <- stringr::str_replace_all(string, "#", "_hash_")
  string <- stringr::str_replace_all(string, "@", "_at_")
  string <- stringr::str_replace_all(string, "\\$", "_dollar_")
  string <- stringr::str_replace_all(string, "\\^", "_caret_")
  string <- stringr::str_replace_all(string, "<", "_less_than_")
  string <- stringr::str_replace_all(string, ">", "_greater_than_")
  string <- stringr::str_replace_all(string, "\\[", "_lbracket_")
  string <- stringr::str_replace_all(string, "\\]", "_rbracket_")
  string <- stringr::str_replace_all(string, "\\{", "_lbrace_")
  string <- stringr::str_replace_all(string, "\\}", "_rbrace_")
  string <- stringr::str_replace_all(string, "\\|", "_pipe_")
  string <- stringr::str_replace_all(string, "\\\\", "_backslash_")
  string <- stringr::str_replace_all(string, "/", "_slash_")
  string <- stringr::str_replace_all(string, "__", "_")  # Remove double underscores
  return(string)
}




#' Map Conditions in Seurat Object
#'
#' This function maps the conditions in the metadata of a Seurat object according to a provided mapping.
#'
#' @param seuratObj A Seurat object.
#' @param metadataName The name of the metadata column in the Seurat object to map from.
#' @param caseCondition The condition value to be mapped to "case".
#' @param controlCondition The condition value to be mapped to "control".
#'
#' @return The Seurat object with the mapped conditions added to its metadata.
#' @importFrom dplyr mutate left_join select
#' @importFrom tibble tibble
#' @export
mapConditionsInSeurat <- function(seuratObj, metadataName, caseCondition, controlCondition) {
  conditionMap <- tibble(
    metadata_values = c(caseCondition, controlCondition),
    condition_mapping = c("case", "control")
  )

  seuratObj@meta.data <- seuratObj@meta.data %>%
    dplyr::mutate(
      condition = dplyr::left_join(
        tibble(metadata_value = .[[metadataName]]),
        conditionMap,
        by = c("metadata_value" = "metadata_values")
      )$condition_mapping
    )

  return(seuratObj)
}

#' Rename Column in a List of Data Frames
#'
#' This function renames a specified column in each data frame within a list of data frames.
#'
#' @param decipher_scores_by_regulon A list where each element is a data frame.
#' @param original_name The original name of the column to be renamed.
#' @param new_name The new name for the column.
#'
#' @return A list of data frames with the specified column renamed in each data frame.
#' @importFrom dplyr rename all_of
#' @export
#' @examples
#' # Assuming 'listDFs' is a list of data frames and you want to rename column 'oldCol' to 'newCol':
#' renamedList <- listOfDFsRenameColumn(listDFs, "oldCol", "newCol")
listOfDFsRenameColumn <- function(decipher_scores_by_regulon,original_name,new_name){
  lookup <- c(original_name)
  names(lookup) <- new_name
  decipher_scores_by_regulon <- decipher_scores_by_regulon %>%
    rename(all_of(lookup))
  return(decipher_scores_by_regulon)
}

#' Add List Element Names to Data Frames Within a List
#'
#' Iterates through a list of data frames and adds a new column to each data frame
#' with the name of the list element.
#'
#' @param this_list A list of data frames.
#' @param element_name The name of the new column to be added to each data frame,
#'        which will contain the name of the list element.
#'
#' @return The same list with each data frame now containing an additional column
#'         named as specified by `elementName`, filled with the name of the list element.
#'
#' @export
#' @examples
#' # Assuming 'listOfDFs' is a list of data frames:
#' updatedList <- addListNameToDFElements(listOfDFs, "sourceName")
addListNameToDFElements <- function(this_list,element_name){
  for(this_element in names(this_list)){
    this_list[[this_element]][[element_name]] <- this_element
  }
  return(this_list)
}



#' Read Dataset parameters
#'
#' @param filePath
#'
#' @return
#' @export
#'
#' @examples
readParameters <- function(filePath) {
  # Check if the file exists
  if (!file.exists(filePath)) {
    stop("File not found: ", filePath)
  }

  # Read the lines from the file
  lines <- readLines(filePath)

  # Initialize an empty list to store parameters
  parameters <- list()

  # Loop through each line and split key and value
  for (line in lines) {
    # Split the line into key and value
    key_value <- strsplit(line, "=")[[1]]
    key <- trimws(key_value[1])
    value <- trimws(key_value[2])

    # Add to the parameters list
    parameters[[key]] <- value
  }

  return(parameters)
}



