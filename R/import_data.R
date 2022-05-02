#' read_xerxes_ras_blocks
#'
#' @param fn path. The path to the input blockTableFile.
#' @param strip_selection_language logical. Should trident selection language marks be stripped from Lefts and Rights?
#' If true, selection marks such as "<" and ">" will be stripped from the Left and Right columns.
#'
#' @return A tibble containing the contents of the blockTableFile
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
read_xerxes_ras_blocks <- function(fn, strip_selection_language = F) {
  if (!file.exists(fn)) {
    stop(paste0(
      "[read_xerxes_ras_blocks] File '",fn, "' not found. Please provide the path to a valid xerxes ras blockTableFile."
    ))
  }

  block_table <- readr::read_tsv(fn, col_types = 'ccicicidd')

  if (strip_selection_language) {
    block_table <- block_table %>%
      dplyr::mutate(dplyr::across(
        .cols=c(.data$Left, .data$Right),
        .fns= ~stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return (block_table)
}
