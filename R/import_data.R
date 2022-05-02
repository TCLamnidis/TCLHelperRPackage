read_xerxes_ras_blocks <- function(fn, stip_selection_language = F) {
  if (!file.exists(fn)) {
    stop(paste0(
      "[read_xerxes_ras_blocks] File '",fn, "' not found. Please provide the path to a valid xerxes ras blockTableFile."
    ))
  }

  block_table <- readr::read_tsv(fn, col_types = 'ccicicidd')

  if (stip_selection_language) {
    block_table <- block_table %>%
      mutate(across(
        .cols=c(Left, Right),
        .fns= ~stringr::str_replace_all(., "([<>])", "") ## Remove individual selection markers
      ))
  }
  return (block_table)
}
