#' Simulate catalogues from multiple signatures
#'
#' This function perfectly combines multiple signatures in proportions described by the `model` argument,
#' then simulates catalogues based on multinomial sampling from the combined signature.
#'
#' Note order of channels in resulting catalogue may differ from their order in
#' the signature collection due to sorting done by [sigstats::sig_combine()].
#'
#' @param signatures A list of signatures to combine and simulate from.
#'   Each signature should be sigverse signature format ([sigshared::example_signature_collection()]).
#' @param model A model for combining the signatures. This can represent the
#'   proportions or contributions of each signature in the combined result.
#' @inheritParams sig_simulate_catalogues_from_signature
#' @inheritParams sigstats::sig_combine
#' @inherit sig_simulate_catalogues_from_signatures return
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load('COSMIC_v3.4_SBS_GRCh38')
#' model <- c(SBS1 = 0.6, SBS5 = 0.4)
#' sig_simulate_catalogues_from_signatures(signatures, model, n = 400, n_catalogues=3)
#'
#' @export
sig_simulate_catalogues_from_signatures <- function(
    signatures, model, n, n_catalogues, format = c("sigverse", "matrix"), prefix = "catalogue", seed = NULL, simplify = FALSE
) {
  # Ensure the 'format' parameter is one of the accepted values
  format <- rlang::arg_match(format)

  # Combine multiple signatures using the provided model
  # This produces a combined signature with contributions specified by 'model'
  combined_sig <- sigstats::sig_combine(signatures = signatures, model = model, format = "signature")


  # Simulate catalogues from the combined signature using the sig_simulate_catalogues_from_signature() function
  sig_simulate_catalogues_from_signature(
    signature = combined_sig,
    n = n,
    n_catalogues = n_catalogues,
    format = format,
    seed = seed,
    prefix = prefix,
    simplify = simplify
  )
}


#' Simulate catalogues from a signature
#'
#' This function performs multinomial sampling from a given mutational signature
#' to generate one or more mutation catalogues. The catalogues can be returned
#' in either 'sigverse' format or as a matrix.
#'
#' @param signature A signature object from the sigverse package. The signature
#'   should contain `fraction` and `channel` columns. See [sigshared::example_signature()].
#' @param n The number of mutations to sample in each catalogue.
#' @param n_catalogues The number of catalogues to generate (default is 1).
#' @param format The format to return the generated catalogues in. Choose between:
#'   - "sigverse": Return a sigverse-style catalogue collection. See [sigshared::example_catalogue_collection()].
#'   - "matrix": Return a matrix where each column corresponds to a simulated catalogue,
#'     with rownames as channels and values as mutation counts.
#' @param seed random seed to use for generation.
#' @param simplify A logical value indicating whether to return a single catalogue
#'   (as a tibble) when only one catalogue is generated and `format = "sigverse"`.
#'   Defaults to FALSE.
#' @param prefix A string prefix to use for naming each catalogue. Defaults to "catalogue".
#'
#' @return Depending on the format:
#'   - If `format = "sigverse"`, a sigverse-style catalogue collection
#'     ([sigshared::example_catalogue_collection()]) is returned.
#'   - If `format = "matrix"`, a matrix with one column per simulated catalogue
#'     and rownames as mutation channels is returned.
#'   - If `simplify = TRUE` and only one catalogue is generated, a single catalogue
#'     tibble ([sigshared::example_catalogue()]) is returned.
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load('COSMIC_v3.4_SBS_GRCh38')
#' sig_simulate_catalogues_from_signature(signatures[["SBS1"]], n = 400, n_catalogues = 3)
#'
#' @export
sig_simulate_catalogues_from_signature <- function(
    signature, n, n_catalogues, format = c("sigverse", "matrix"),
    prefix = "catalogue", seed = NULL, simplify = FALSE
) {
  # Ensure the 'format' parameter is one of the accepted values
  format <- rlang::arg_match(format)

  # Assert that the provided signature is valid
  sigshared::assert_signature(signature)

  # Extract fractions and channels from the signature
  fractions <- signature[["fraction"]]
  names(fractions) <- signature[["channel"]]


  # Simulate catalogues using multinomial distribution
  # `rmultinom` generates counts for each channel across multiple catalogues
  mx <- sigshared::with_seed(
    seed = seed,
    {rmultinom(n = n_catalogues, size = n, prob = fractions)}
  )

  # Assign names to each catalogue using the provided prefix
  colnames(mx) <- paste0(prefix, "_", seq_len(ncol(mx)))

  # Return as matrix if format is "matrix"
  if (format == "matrix") {
    return(mx)
  }

  # Setup the catalogue structure with 'channel' and 'type' columns
  catalogue_std_columns <- sigshared::bselect(signature, c("channel", "type"))

  # Split the simulated catalogues (matrix) into a list (by columns) for easier processing
  ls_counts <- asplit(mx, MARGIN = 2)

  # Create sigverse-style catalogues (one element for each simulated catalogue)
  catalogues <- lapply(ls_counts, function(count) {
    catalogue <- cbind(catalogue_std_columns, count = count)

    # Calculate fraction of each mutation type in the catalogue
    catalogue[["fraction"]] <- catalogue[["count"]] / sum(catalogue[["count"]])

    tibble::tibble(catalogue)
  })

  # If only one catalogue is generated and simplify = TRUE, return a single catalogue
  if (simplify && length(ls_counts) == 1) {
    return(catalogues[[1]])
  }

  # Otherwise, return a collection of catalogues
  return(catalogues)
}

perfect_reconstruction <- function(signature, n, n_catalogues=1){
  n_channels = nrow(signature)
  mx = matrix(rep(signature[["fraction"]] * n , times = n_catalogues), nrow = n_channels)
}
