test_that("I can produce a single meta-cell", {
  #input data
  rnaCountsMatrix <- matrix(c(1,2,1,3,1,2,1,1,4),nrow=3)
  colnames(rnaCountsMatrix) <- c("A","B","C")
  distanceMatrix <- CorDist(rnaCountsMatrix)
  numNearestNeighbors <- 2
  numMetaCells <- 1

  decipher_result <- calculatePseudoBulkMatrix(rnaCountsMatrix, distanceMatrix, numNearestNeighbors, numMetaCells)
  decipher_result <- c(decipher_result)
  expected_result <- rowSums(rnaCountsMatrix)
  expect_equal(decipher_result, expected_result)
})

test_that("function does not include any one cell more than once", {
  #input data
  rnaCountsMatrix <- matrix(c(1,2,1,3,1,2,1,1,4),nrow=3)
  colnames(rnaCountsMatrix) <- c("A","B","C")
  distanceMatrix <- CorDist(rnaCountsMatrix)
  numNearestNeighbors <- 3
  numMetaCells <- 1
  decipher_result <- calculatePseudoBulkMatrix(rnaCountsMatrix, distanceMatrix, numNearestNeighbors, numMetaCells)
  decipher_result <- c(decipher_result)
  expected_result <- rowSums(rnaCountsMatrix)
  expect_equal(decipher_result, expected_result)
})

test_that("function does not return more meta cells than possible (no overlap)", {
  #input data
  rnaCountsMatrix <- matrix(c(1,2,1,3,1,2,1,1,4),nrow=3)
  colnames(rnaCountsMatrix) <- c("A","B","C")
  distanceMatrix <- CorDist(rnaCountsMatrix)
  numNearestNeighbors <- 2
  numMetaCells <- 2
  decipher_result <- calculatePseudoBulkMatrix(rnaCountsMatrix, distanceMatrix, numNearestNeighbors, numMetaCells)
  decipher_result <- c(decipher_result)
  expected_result <- "warning, more meta-cells requested than possible without overlap"
  expect_equal(decipher_result, expected_result)
})

test_that("don't return overlapping meta cells", {
  #input data
  rnaCountsMatrix <-matrix(sample(1:100, 100*100, replace = TRUE), nrow = 100, ncol = 100)

  # Create column and row names
  col_names <- paste("cell", 1:100)
  row_names <- paste("gene", 1:100)

  # Set the column and row names for the matrix
  colnames(rnaCountsMatrix) <- col_names
  rownames(rnaCountsMatrix) <- row_names

  seurat_object <- Seurat::CreateSeuratObject(counts = rnaCountsMatrix)
  seurat_object@meta.data$cluster <- "this_cluster"
  seurat_object@meta.data$condition <- rep(c("case","control"),each = 50)
  decipher_result <- generateMetaCellMatrices(seurat_object,paramMinMetaCells = 1, paramMaxMetaCells = 600, paramMaxScCells = 3000, paramK = 2)
  actual_result <- dim(decipher_result$this_cluster$case)[2]
  expected_result <- 25
  expect_equal(actual_result, expected_result)
})




