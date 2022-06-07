#' ras_to_similarity_matrix
#'
#' @param ras_table tibble. A tibble containing the contents of the tableOutFile.
#' @param na_diag logical. If TRUE the diagonal of the similarity matrix is set to NA.
#'
#' @return tibble. A similarity matrix where the Right populations are columns and their values are the ras of each Left with each right. Standard errors are discarded.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data

ras_to_similarity_matrix <- function(ras_table, na_diag=F) {
  similarity_matrix <- ras_table %>%
    tidyr::pivot_wider(id_cols=.data$Left, names_from = .data$Right, values_from = .data$RAS)

  row_names <- similarity_matrix$Left
  similarity_matrix <- similarity_matrix %>%
    dplyr:: select(-.data$Left) %>%
    as.matrix()
  rownames(similarity_matrix) <- row_names

  if ( na_diag ) {
    diag(similarity_matrix) <- NA
  }

  return(similarity_matrix)
}
