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
      "[read_xerxes_ras_blocks] File '", fn, "' not found. Please provide the path to a valid xerxes ras blockTableFile."
    ))
  }

  block_table <- readr::read_tsv(fn, col_types = "ccicicidd")

  if (strip_selection_language) {
    block_table <- block_table %>%
      dplyr::mutate(dplyr::across(
        .cols = c(.data$Left, .data$Right),
        .fns = ~ stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return(block_table)
}

#' read_xerxes_ras_table
#'
#' @param fn path. The path to the input tableOutFile.
#' @param strip_selection_language logical. Should trident selection language marks be stripped from Lefts and Rights?
#' If true, selection marks such as "<" and ">" will be stripped from the Left and Right columns.
#'
#' @return A tibble containing the contents of the tableOutFile.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
read_xerxes_ras_table <- function(fn, strip_selection_language = F) {
  if (!file.exists(fn)) {
    stop(paste0(
      "[read_xerxes_ras_table] File '", fn, "' not found. Please provide the path to a valid xerxes ras tableOutFile"
    ))
  }

  ras_table <- readr::read_tsv(fn, col_types = "ccddd")

  if (strip_selection_language) {
    ras_table <- ras_table %>%
      dplyr::mutate(dplyr::across(
        .cols = c(.data$Left, .data$Right),
        .fns = ~ stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return(ras_table)
}

#' read_xerxes_f3_table
#'
#' @param fn path. The path to the input tableOutFile.
#' @param strip_selection_language logical. Should trident selection language marks be stripped from columns a,b and c?
#' If true, selection marks such as "<" and ">" will be stripped from the a, b and c columns.
#'
#' @return A tibble containing the contents of the tableOutFile.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
read_xerxes_f3_table <- function(fn, strip_selection_language = F) {
  if (!file.exists(fn)) {
    stop(paste0(
      "[read_xerxes_f3_table] File '", fn, "' not found. Please provide the path to a valid xerxes f3 tableOutFile"
    ))
  }
  ## Skip header and add own column names because of https://github.com/poseidon-framework/poseidon-analysis-hs/issues/3
  f3_table <- readr::read_tsv(fn, col_types = "cccccddd", skip=1, col_names=c("Statistic", "a", "b", "c", "Estimate","StdErr","Z score"))

  if (strip_selection_language) {
    f3_table <- f3_table %>%
      dplyr::mutate(dplyr::across(
        .cols = c(.data$a, .data$b, .data$c),
        .fns = ~ stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return(f3_table)

}
