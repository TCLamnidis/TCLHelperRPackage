#' f3_to_similarity_matrix
#'
#' @param f3_table tibble. A tibble containing the contents of the F3tableOutFile.
#' @param na_diag logical. If TRUE the diagonal of the similarity matrix is set to NA.
#'
#' @return tibble. A similarity matrix where the 'b' populations are columns and their values are the F3 of each 'a' with each 'b'. Standard errors are discarded.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data

f3_to_similarity_matrix <- function(f3_table, na_diag=F) {
  ## First, check that the dataset only contains one outgroup, else a similarity matrix makes no sense.
  if ( f3_table %>% dplyr::count(c) %>% nrow() != 1 ) {
    stop(paste0(
      "[f3_to_similarity_matrix] Cannot create a similarity matrix with statistics using multiple outgroups."
    ))
  }

  similarity_matrix <- f3_table %>%
    tidyr::pivot_wider(id_cols=.data$a, names_from = .data$b, values_from = .data$Estimate)

  row_names <- similarity_matrix$a
  similarity_matrix <- similarity_matrix %>%
    dplyr:: select(-.data$a) %>%
    as.matrix()
  rownames(similarity_matrix) <- row_names

  if ( na_diag ) {
    diag(similarity_matrix) <- NA
  }

  return(similarity_matrix)
}
