#' Delete mj jackknife function
#'
#' Implementation of delete mj jackknifing following Nick Patterson's writeup of it.
#' The function works like tidyverse functions, i.e. values, mj and weights should not be quoted/be symbols.
#'
#' @param data Tibble containing the data to use.
#' @param values Unquoted name of the tibble column with the block values.
#' @param mj Unquoted name of the tibble column with the mj values.
#' @param weights Unquoted name of the tibble column with the block weights.
#' @param .print_vals logical. Should the actual theta_minus values be printed?
#'
#' @return list. Returns a list of items :
#' theta_J  The jackknife mean from all the theta_-x data.
#' jack_se The SE of theta_J
#' Zscore  The Zscore of the theta_J
#' jack.values A vector of all the theta_-x values
#' @export
#'
delete_mj_jackknife <- function(data, values, mj, weights, .print_vals=F) {
  values <- rlang::enquo(values)
  mj <- rlang::enquo(mj)
  weights <- rlang::enquo(weights)

  data <- data %>%  dplyr::ungroup() ## Remove grouping from the data to avoid summarise mishaps

  ## Calculate the sum of each column as well as theta_hat
  sums <- data %>%  dplyr::summarise(value_sum=sum(!!values),     ## The sum of values
                                     mj_sum=sum(!!mj),            ## The sum of mj. used for calculating rates of the values
                                     weight_sum=sum(!!weights),   ## The sum of weights across all blocks
                                     theta_hat=.data$value_sum/.data$mj_sum   ## Theta hat as a rate of values across all blocks
  )

  ## Add columns with the theta_minus, as well as n_minus_weight
  data <- data %>% dplyr::mutate(
    value_minus=sums$value_sum-!!values,        ## value total - value of block
    mj_minus=sums$mj_sum-!!mj,                  ## mj total - mj of block
    n_minus_weight=sums$weight_sum-!!weights,   ## weight total - weight of block
    theta_minus=.data$value_minus/.data$mj_minus,           ## The rate calculated from all but the given block.
    hj=sums$weight_sum/!!weights,               ## hj per block for tau calculation
    tau=.data$hj*sums$theta_hat - (.data$hj-1)* .data$theta_minus ## pseudovalue for variance calculation
  )
  ## Use added columns to calculate theta_J
  theta_J <- data %>%
    dplyr::summarise(
      theta_J=nrow(data) * sums$theta_hat - sum(.data$n_minus_weight*.data$theta_minus/sums$weight_sum)
      ) %>%
    dplyr::pull(.data$theta_J)
  variance <- data %>%
    dplyr::summarise(
      variance=1/nrow(data) * sum((.data$tau-theta_J)^2/(.data$hj-1))
      ) %>%
    dplyr::pull(.data$variance)
  jack_se = sqrt(variance)
  if (.print_vals) {
    return(
      list(
        theta_J=theta_J,
        jack_se=jack_se,
        Zscore=theta_J/jack_se,
        jack.values=as.vector(data$theta_minus)
      )
    )
  } else {
    return(
      list(
        theta_J=theta_J,
        jack_se=jack_se,
        Zscore=theta_J/jack_se
      )
    )
  }
}
