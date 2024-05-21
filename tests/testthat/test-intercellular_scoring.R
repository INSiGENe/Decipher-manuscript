test_that("multiplication works", {
  receptorMatrix <-   matrix(sample(1:100, 2*2, replace = TRUE), nrow = 2, ncol = 2)

  # Create column and row names
  col_names <- paste("cell", 1:2)
  row_names <- paste("gene", 1:2)

  # Set the column and row names for the matrix
  colnames(receptorMatrix) <- col_names
  rownames(receptorMatrix) <- row_names

  conditionVector <- c("case","control")
  names(conditionVector) <- paste("cell", 1:2)

  ligandMeans <- data.frame(
    ligand = paste("gene", 3:4),
    case = c(0,1),
    control = c(1,0)
  )
  rownames(ligandMeans) <- ligandMeans$ligand
  LRSet <- data.frame(
    receptor = c("gene 1","gene 2"),
    ligand = c("gene 3","gene 4"),
    interaction = c("gene 3-gene 1","gene 4-gene 2")
  )

  decipher_result <- calculateInteractionMatrix(receptorMatrix, conditionVector, ligandMeans, LRSet)
  decipher_result <- decipher_result^2
  expected_result <- matrix(
    0,
    nrow = 2, ncol = 2)

  expected_result [1,2] <- receptorMatrix[1,2]
  expected_result [2,1] <- receptorMatrix[2,1]
  colnames(expected_result) <- c("cell 1","cell 2")
  rownames(expected_result) <- c("gene 3-gene 1","gene 4-gene 2")
  expect_equal(decipher_result, expected_result)
})
