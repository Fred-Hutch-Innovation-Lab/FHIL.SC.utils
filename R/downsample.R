#' Downsample counts data
#'
#' Downsample counts matrices to target read counts using multinomial sampling.
#' Downsampling can be performed at the whole matrix level, or at the
#' column/cell level by setting `per_cell = TRUE` (the default). If the input
#' is a Seurat object, the function will return a Seurat object with some optional
#' tidying if `tidy_seurat_obj = TRUE`. The tidying may cause issues if Seurat
#' object organization changes, so you may want disable this if you get errors.
#'
#' @param data A (sparse) counts matrix or Seurat Object
#' @param nreads Target number of reads to downsample to
#' @param per_cell Is the target reads per cell, or for the whole dataset?
#' @param tidy_seurat_obj Flag to do additional calculations and tidying of
#'  Seurat output based on new values
#'
#' @import Matrix
#' @importFrom stats rmultinom
#'
#' @return A Seurat object or matrix, echoing the class of the input
#' @export
#'
#' @examples {
#'
#' data <- FHIL.SC.utils::pbmc3k
#' downsampled_data <- downsample(data, nreads=10000, tidy_seurat_obj=TRUE)
#'
#' head(data$nCount_RNA)
#' head(downsampled_data$nCount_RNA)
#'
#' }
downsample <- function(data, nreads, per_cell=TRUE, tidy_seurat_obj=TRUE) {
  if (inherits(data, 'Seurat')) {
    counts <- data@assays$RNA@layers$counts
  } else if (inherits(data, 'matrix') | inherits(data, 'dgCMatrix')) {
    counts <- data
  } else {
    stop('Data input should be a matrix, sparse matrix, or Seurat Object')
  }

  downsample_vector <- function(count_vector, target_sum) {
    original_sum <- sum(count_vector)
    probs <- count_vector / original_sum
    downsampled_counts <- rmultinom(1, size = target_sum, prob = probs)
    return(as.vector(downsampled_counts))
  }

  if (per_cell) {
    downsampled_data <- apply(counts, 2, function(col) downsample_vector(col, nreads))
  } else {
    downsampled_data <- downsample_vector(counts, nreads)
  }

  ## This seems like overkill, but it will be preserved here in case it's wanted
  # # Make downsampling exactly to nreads
  # diff <- nreads - sum(downsampled_data)
  #
  # if (diff != 0) {
  #   # Find the non-zero entries that can be adjusted
  #   non_zero_indices <- which(as.vector(counts > 0))
  #
  #   # Use the original probabilities to weight the adjustments
  #   adjust_probs <- probs[non_zero_indices]
  #   adjust_probs <- probabilities[non_zero_indices]
  #
  #   # Normalize the adjustment probabilities to sum to 1
  #   adjust_probs <- adjust_probs / sum(adjust_probs)
  #
  #   # Randomly select indices to adjust based on the adjusted probabilities
  #   if (diff > 0) {
  #     # If we need to add counts, increment random indices weighted by original probabilities
  #     adjust_indices <- sample(non_zero_indices, diff, replace = TRUE, prob = adjust_probs)
  #     downsampled_counts[adjust_indices] <- downsampled_counts[adjust_indices] + 1
  #   } else {
  #     # If we need to reduce counts, decrement random non-zero indices weighted by original probabilities
  #     adjust_indices <- sample(non_zero_indices, abs(diff), replace = TRUE, prob = adjust_probs)
  #     downsampled_counts[adjust_indices] <- pmax(downsampled_counts[adjust_indices] - 1, 0)
  #   }
  # }

  # reform matrix from vector
  downsampled_matrix <- Matrix::Matrix(data=matrix(downsampled_data,
                                                   nrow=nrow(counts),
                                                   ncol=ncol(counts)),
                                       sparse=TRUE)

  if (all(dim(counts) != dim(downsampled_matrix))) {
    stop('Error creating downsampled matrix: output is not of same dimensions as input')
  }

  if (inherits(data, 'Seurat')) {
    if (tidy_seurat_obj) {
      if (length(names(data@assays$RNA@layers)) > 1) {
        message('Removing additional layers besides counts.
              Normalized and/or scaled matrices will have to be recalculated')
      }
      data@assays$RNA@layers <- list()
      data@assays$RNA@layers$counts <- downsampled_matrix

      data$nCount_RNA <- colSums(data@assays$RNA@layers$counts)
      data$nFeature_RNA <- colSums(data@assays$RNA@layers$counts > 0)


    }
  } else if (inherits(data, 'matrix') | inherits(data, 'dgCMatrix')) {
    data <- downsampled_matrix
  }
  return(data)
}
