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
#' @inheritParams read_xerxes_ras_blocks
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

#' read_xerxes_fstats_table
#'
#' Read xerxes tableOutFile into a tibble. F2 statistics are converted to F3 statistics to simplify the creation of a similarity matrix.
#'
#' @param fn path. The path to the input tableOutFile.
#' @inheritParams read_xerxes_ras_blocks
#' @param old_xerxes logical. Older versions of xerxes did not correctly match the columns for F3/F2 to the header (which is set up for F4).
#' Should be set to FALSE for tables produced with xerxes version higher than v0.1.3.1.
#'
#' @return A tibble containing the contents of the tableOutFile.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
read_xerxes_fstats_table <- function(fn, strip_selection_language = F, old_xerxes=F) {
  if (!file.exists(fn)) {
    stop(paste0(
      "[read_xerxes_f3_table] File '", fn, "' not found. Please provide the path to a valid xerxes f3 tableOutFile"
    ))
  }

  if (old_xerxes) {
    ## Skip header and add own column names because of https://github.com/poseidon-framework/poseidon-analysis-hs/issues/3
    f3_table <- data.table::fread(fn, sep="\t", header=T, fill=T) %>%
      homogenise_xerxes_fstats_output(., old_xerxes)
  } else {
    f3_table <- data.table::fread(fn, sep="\t", header=T, fill=T) %>%
      dplyr::mutate(c=dplyr::na_if(.data$c, '')) %>%
      homogenise_xerxes_fstats_output(., old_xerxes)
  }

  if (strip_selection_language) {
    f3_table <- f3_table %>%
      dplyr::mutate(dplyr::across(
        .cols = c(.data$a, .data$b, .data$c),
        .fns = ~ stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return(f3_table)

}

homogenise_xerxes_fstats_output <- function(input_f3_table, old_xerxes) {
    f3_table <- input_f3_table %>%
      dplyr::mutate(n_na=rowSums(is.na(.)))
    if (old_xerxes) {
      f2_table <- f3_table %>%
        dplyr::filter(.data$n_na == 2) %>%
        dplyr::mutate(
          `Z score` = .data$Estimate,
          StdErr = as.double(.data$d),
          Estimate = as.double(.data$c),
          d = NA_character_,
          c = .data$b,
          b = .data$a
          )
      f4_table <- f3_table %>%
        dplyr::filter(.data$n_na == 0) %>%
        dplyr::mutate(
          dplyr::across(
            .cols=c("Statistic", "a", "b", "c", "d"), .fn=~as.character()
          ),
          dplyr::across(
            .col=c("Estimate", "StdErr", "Z score"), .fn=~as.numeric()
          )
        )
      f3_table <- f3_table %>%
        dplyr::filter(.data$n_na == 1) %>%
        dplyr::mutate(
          `Z score` = .data$StdErr,
          StdErr = .data$Estimate,
          Estimate = as.double(.data$d),
          d = NA_character_
        )
      homogenised_table <- dplyr::bind_rows(f2_table, f3_table, f4_table) %>%
        dplyr::select(-.data$n_na) %>%
        tibble::as_tibble()
    } else {
      ## With newer versions of xerxes only f2 needs to be adjusted to F3. F3 results are read correctly.
      f2_table <- f3_table %>%
        dplyr::filter(.data$n_na == 2) %>%
        dplyr::mutate(
          d = NA_character_,
          c = .data$b,
          b = .data$a
        )
      f3_table <- f3_table %>%
        dplyr::filter(! .data$n_na == 2)

      homogenised_table <- dplyr::bind_rows(f2_table, f3_table) %>%
        dplyr::select(-.data$n_na) %>%
        tibble::as_tibble()
    }

    return(homogenised_table)
}

#' read_xerxes_rasf4_table
#'
#' Read xerxes ras f4TableOutFile contents into a tibble.
#'
#' @param fn character. The path to the input f4TableOutFile.
#' @inheritParams read_xerxes_ras_blocks
#'
#' @return A tibble containing the contents of the f4TableOutFile.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
read_xerxes_rasf4_table <- function(fn, strip_selection_language = F) {
  rasf4_table <- readr::read_tsv(fn, col_types='cccddddd')

  if ( strip_selection_language ) {
    rasf4_table <- rasf4_table %>%
      dplyr::mutate(dplyr::across(
        .cols = c(.data$Left1, .data$Left2, .data$Right),
        .fns = ~ stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }

  return(rasf4_table)
}
