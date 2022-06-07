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
  f3_table <- data.table::fread(fn, sep="\t", header=T, fill=T) %>%
    homogenise_xerxes_output()

  if (strip_selection_language) {
    f3_table <- f3_table %>%
      dplyr::mutate(dplyr::across(
        .cols = c(.data$a, .data$b, .data$c),
        .fns = ~ stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return(f3_table)

}

homogenise_xerxes_output <- function(input_f3_table) {
    f3_table <- input_f3_table %>%
      dplyr::mutate(n_na=rowSums(is.na(.)))
    f2_table <- f3_table %>%
      dplyr::filter(n_na == 2) %>%
      dplyr::mutate(
        `Z score` = .data$Estimate,
        StdErr = as.double(.data$d),
        Estimate = as.double(.data$c),
        d = NA_character_,
        c = .data$b,
        b = .data$a
        )
    f3_table <- f3_table %>%
      dplyr::filter(n_na == 1) %>%
      dplyr::mutate(
        `Z score` = .data$StdErr,
        StdErr = .data$Estimate,
        Estimate = as.double(.data$d),
        d = NA_character_
      )

    homogenised_table <- dplyr::bind_rows(f2_table, f3_table) %>%
      dplyr::select(-n_na) %>%
      tibble::as_tibble()

    return(homogenised_table)
}
