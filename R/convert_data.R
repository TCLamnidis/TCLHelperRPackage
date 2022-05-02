if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Norm","rasta")) ## Disables notes about '.' due to magrittr

#' A function to calculate RASTA(A, B; C, Ref) from a blockTable of RAS per chromosome.
#'
#' @param block_table tibble. A tibble with the contents of a blockTableFile, derived with read_xerxes_ras_blocks()
#' @param A character. The Lefts to cycle as population A.
#' @param B character. The Lefts to cycle as population B.
#' @param C character. The Rights to cycle as population C.
#'
#' @return A tibble with the all the requested combinations of rasta, along with their SE and Zscores
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
ras_to_rasta <- function(block_table, A, B, C) {

  ## Check that all populations were provided, and requested populations exist in the input
  check_x_in_y(A, block_table$Left, "A")
  check_x_in_y(B, block_table$Left, "B")
  check_x_in_y(C, block_table$Right, "C")

  ## Convert ras to raw ras
  raw_ras <- block_table %>%
    dplyr::mutate(RAS=.data$RAS*.data$Norm) %>%
    dplyr::select(.data$Left, .data$Right, .data$BlockNr, .data$RAS, .data$Norm)

  tidyr::crossing(A,B,C) %>% purrr::pmap_dfr(., ~calculate_rasta(raw_ras, A=..1, B=..2, C=..3))
}


#' Calculate a specific RASTA(A, B, C, Ref) from an input block table.
#'
#' @param raw_ras tibble. A tibble containing the results
#' @param A character. Population A fr the rasta calculation. Must be present in the raw_ras$Left.
#' @param B character. Population B for the rasta calculation. Must be present in the raw_ras$Left.
#' @param C character. Population C for the rasta calculation. Must be present in the raw_ras$Right.
#'
#' @return A tibble with the requested RASTA, it's error and ZScore
#' @export
calculate_rasta <- function(raw_ras, A, B, C) {

  ras_A_C <- raw_ras %>%
    dplyr::filter(.data$Left == A, .data$Right == C)

  ras_B_C <- raw_ras %>%
    dplyr::filter(.data$Left == B, .data$Right == C)

  rasta_table <- dplyr::full_join(ras_A_C, ras_B_C, by=c("Right", "BlockNr","Norm")) %>%
    dplyr::rename(A=.data$Left.x, B=.data$Left.y, C=.data$Right) %>%
    dplyr::mutate(rasta=.data$RAS.x - .data$RAS.y, D="Ref") %>%
    dplyr::select(.data$A, .data$B, .data$C, .data$D, .data$BlockNr, .data$Norm, .data$rasta)

  result <- delete_mj_jackknife(rasta_table, rasta, Norm, Norm) %>% tibble::as_tibble() %>%
    dplyr::mutate(A=A, B=B, C=C, D="Ref") %>%
    dplyr::select(.data$A, .data$B, .data$C, .data$D, .data$theta_J, .data$jack_se, .data$Zscore)

  return(result)
}

#' Check that x is not empty and all elements exist in y.
#'
#' @param x character. A vector of strings that must exist in y.
#' @param y character. A vector of character that must include all elements in x.
#' @param var_name character. The variable name to be included in error messages.
#'
#' @export
#'
check_x_in_y <- function(x, y, var_name) {
  assertthat::assert_that(length(x) != 0, msg = paste0("[ras_to_rasta]: '", var_name, "' cannot be empty."))

  assertthat::assert_that(
    setequal(x, intersect(x, y))
    , msg = paste0("[ras_to_rasta]: ID(s) '", paste0(setdiff(x, y), collapse="', "), "' not found in input table. Halting execution.")
  )
}
