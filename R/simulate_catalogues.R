#' Simulate catalogues from multiple signatures
#'
#' Simulates mutation catalogues by first combining multiple
#' signatures (according to specified proportions in the `model`), and then sampling mutations
#' from the resulting *combined* signature using a multinomial distribution.
#'
#' This approach introduces **realistic sampling noise**, which mimics how mutations
#' arise in actual biological samples. It is ideal for benchmarking how well tools
#' can detect signatures when contributions are probabilistic rather than exact.
#'
#' Unlike [`sig_simulate_stratified()`], this method does **not** guarantee exact
#' counts per signature in the output — the breakdown is subject to sampling variation.
#' To simulate catalogues with *deterministic* signature contributions (no noise),
#' use [`sig_simulate_perfect()`].
#'
#' @param signatures A list of signatures to combine and simulate from.
#'   Each signature should be sigverse signature format ([sigshared::example_signature_collection()]).
#' @param model A model for combining the signatures. This can represent the
#'   proportions or contributions of each signature in the combined result.
#' @inheritParams sig_simulate_single
#' @inheritParams sigstats::sig_combine
#' @inherit sig_simulate_mixed return
#'
#' @seealso
#' - [sig_simulate_perfect()] for noiseless catalogue generation
#' - [sig_simulate_stratified()] for sampling each signature independently
#' - [sig_simulate_single()] for sampling from a single signature only
#'
#' @examples
#' library(sigstash)
#'
#' # Load Signature Collection
#' signatures <- sig_load('COSMIC_v3.4_SBS_GRCh38')
#'
#' # Define Model
#' model <- c(SBS1 = 0.6, SBS5 = 0.4)
#'
#' # Simulate randomly from a mixed signature
#' sig_simulate_mixed(signatures, model, n = 400, n_catalogues=3)
#'
#' @export
sig_simulate_mixed <- function(
    signatures, model, n, n_catalogues, format = c("sigverse", "matrix"), prefix = "catalogue", seed = NULL, simplify = FALSE
) {

  # Ensure the 'format' parameter is one of the accepted values
  format <- rlang::arg_match(format)

  # Combine multiple signatures using the provided model
  # This produces a combined signature with contributions specified by 'model'
  combined_sig <- sigstats::sig_combine(signatures = signatures, model = model, format = "signature")


  # Simulate catalogues from the combined signature using the sig_simulate_single() function
  sig_simulate_single(
    signature = combined_sig,
    n = n,
    n_catalogues = n_catalogues,
    format = format,
    seed = seed,
    prefix = prefix,
    simplify = simplify
  )
}

#' Simulate Signatures by Stratified Sampling
#'
#' Simulate catalogues by independently sampling from each signature in the model,
#' ensuring exact contribution proportions.
#'
#' @inheritParams sig_simulate_mixed
#' @param verbose verbose mode (flag)
#' @return A named list with:
#'  - `mutations_per_signature`: A named vector of how many mutations were sampled from each signature
#'  - `catalogues`: A list or matrix of simulated catalogues (see `format` argument)
#'
#' @examples
#' library(sigstash)
#'
#' #' # Load Signature Collection
#' signatures <- sig_load('COSMIC_v3.4_SBS_GRCh38')
#'
#' # Define Model
#' model <- c(SBS1 = 0.6, SBS5 = 0.4)
#'
#' # Independently sample 0.6x400 mutations from SBS1 & 0.4*400 mutations from SBS5
#' simulations <- sig_simulate_stratified(signatures, model, n = 400, n_catalogues=3)
#'
#' # Confirm how many mutations were are from each signature
#' simulations$mutations_per_signature
#'
#' simulations$catalogues
#'
#' @export
sig_simulate_stratified <- function(
    signatures, model, n, n_catalogues, format = c("sigverse", "matrix"), prefix = "catalogue", seed = NULL, simplify = FALSE, verbose=FALSE
  ){

  # Ensure the 'format' parameter is one of the accepted values
  format <- rlang::arg_match(format)

  # Compute how many mutations from each signature.
  mutations_per_sig = round(model * n, digits = 0)
  signames = names(model)

  if(verbose){
    string_counts <- toString(paste0(names(mutations_per_sig), ": ", mutations_per_sig, " mutations"))
    cli::cli_alert_info("Simulating [{string_counts}]")
  }

   cataloges_per_sig <- lapply(
    seq_along(mutations_per_sig),
    function(i){
      signame = signames[i]
      n_muts = mutations_per_sig[i]
      sig_simulate_single(
        signature = signatures[[signame]], n = n_muts, n_catalogues = n_catalogues, prefix = signame
      )
    }
  )

   names(cataloges_per_sig) <- signames

   catalogues_summed <- lapply(
     seq_len(n_catalogues),
     function(i){
       cats_to_sum <- lapply(cataloges_per_sig, function(c) { c[[i]] })
       sigstats::sig_sum(cats_to_sum)
     }
    )

   names(catalogues_summed) <- paste0(prefix, "_", seq_along(catalogues_summed))

   if(format == "matrix"){
     catalogues_summed <- sigshared::sig_collection_reformat_list_to_matrix(catalogues_summed, values = "count")
   }

   # Add indication of how many mutations are from each signature
   ls_simulations <- list(
     mutations_per_signature = as.list(mutations_per_sig),
     catalogues = catalogues_summed
   )

   return(ls_simulations)

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
#' sig_simulate_single(signatures[["SBS1"]], n = 400, n_catalogues = 3)
#'
#' @export
sig_simulate_single <- function(
    signature, n, n_catalogues, format = c("sigverse", "matrix"),
    prefix = "catalogue", seed = NULL, simplify = FALSE
) {
  # Ensure the 'format' parameter is one of the accepted values
  format <- rlang::arg_match(format)

  # Assertions
  sigshared::assert_signature(signature)
  assertions::assert_greater_than(n_catalogues, minimum = 0)
  assertions::assert_greater_than_or_equal_to(n, minimum = 0)

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

  # Add type attribute (vector of types matched to channel names)
  attr(mx, "type") <- signature[["type"]]

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

#' Simulate a Perfect (Noiseless) Catalogue from a Signature Model
#'
#' This function creates a single synthetic mutation catalogue by \strong{deterministically} combining
#' multiple mutational signatures in fixed proportions (as specified by a model), then scaling
#' the resulting signature to exactly `n` mutations.
#' No randomness is introduced—this is a "perfect" profile.
#'
#' This is useful for generating ground truth references or benchmarking results without the noise introduced
#' by multinomial sampling (e.g. as used in `sig_simulate_mixed()`).
#'
#' @param signatures A named list of signatures in sigverse format. See [sigshared::example_signature_collection()].
#' @param model A named numeric vector representing the proportion of each signature to include in the combined profile.
#'   The names must match entries in `signatures`. See [sigshared::example_model()]
#' @param n The total number of mutations to simulate.
#' @param format The output format. Choose between:
#'   - `"sigverse"`: Returns a single sigverse-style catalogue data.frame with `channel`, `type`, `count`, and `fraction`.
#'   - `"matrix"`: Returns a matrix with `channel` names as rows and counts as values.
#' @param round round counts to ensure they're whole numbers.
#' This will introduce some 'imperfection' in that fraction of counts may differ from underlying signature.
#' @return A single perfect catalogue representing the deterministic combination of signatures in `model`.
#'   Returned in either sigverse-style format or matrix form, depending on `format` argument.
#'
#' @seealso [sig_simulate_mixed()] for noisy simulation from mixed signatures,
#'   [sig_simulate_stratified()] for stratified sampling, and [sigstats::sig_reconstruct()] for
#'   noiseless single-signature catalogues.
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load("COSMIC_v3.4_SBS_GRCh38")
#' model <- c(SBS1 = 0.7, SBS5 = 0.3)
#'
#' perfect <- sig_simulate_perfect(signatures, model = model, n = 400)
#' head(perfect)
#'
#' # Return as matrix
#' perfect_mat <- sig_simulate_perfect(signatures, model = model, n = 400, format = "matrix")
#'
#' @export
sig_simulate_perfect <- function(signatures, model, n, format = c("sigverse", "matrix"), round=TRUE){

  # Assertions
  sigshared::assert_signature_collection(signatures)
  sigshared::assert_model(model)
  assertions::assert_number(n)
  format <- rlang::arg_match(format)

  # Create combined sig based on model
  sig_combined <- sigstats::sig_combine(signatures = signatures, model = model, format = "signature", verbose = FALSE)

  # Create Counts
  sig_combined[["count"]] <- sig_combined[["fraction"]] * n

  # Round counts
  if(round){
    sig_combined[["count"]] <- round(sig_combined[["count"]], digits=0)

    # Fix fraction of catalogue (can be changed by rounding)
    sig_combined[["fraction"]] <- sigshared::compute_fraction(sig_combined[["count"]])
  }

  # Move fraction to last column
  catalogue <- sig_combined[c(setdiff(colnames(sig_combined), "fraction"), "fraction")]

  # Optionally Reformat to matrix
  if(format == "matrix"){
    catalogue <- sigshared::sig_collection_reformat_list_to_matrix(list("simulated"=catalogue), values = "count")
  }

  return(catalogue)
}



#' Bootstrap a Mutation Catalogue
#'
#' Generates bootstrap resamples of a mutational catalogue by multinomial resampling
#' of its observed mutation counts. This allows creation of synthetic catalogues that
#' reflect sampling variability in mutation observations.
#'
#' The probability distribution for resampling is based on the `fraction` or `count` column
#' of the input catalogue. A total of `n_catalogues` bootstrapped catalogues are created,
#' each with the same total number of mutations as the original catalogue.
#'
#' @param catalogue A single sigverse-style catalogue to bootstrap. See [sigshared::example_catalogue()].
#' @param n_catalogues The number of bootstrap replicates to generate. Default is 10.
#' @param format The output format. Either:
#'   - `"sigverse"`: return a list of sigverse-style data.frames (default).
#'   - `"matrix"`: return a matrix with rows as channels and columns as simulated catalogues.
#' @param prefix String prefix to name each catalogue. Default is `"catalogue"`.
#' @param seed Optional seed for reproducibility.
#'
#' @return Either a list of sigverse-style catalogues or a matrix of bootstrapped counts.
#'
#' @examples
#' library(sigshared)
#'
#' # Bootstrap catalogue (return catalogue collection as list)
#' bootstraps <- sig_bootstrap_catalogue(example_catalogue(), n_catalogues = 5)
#'
#' # Return  matrix
#' mx <- sig_bootstrap_catalogue(example_catalogue(), n_catalogues = 5, format = "matrix")
#'
#' @export
sig_bootstrap_catalogue <- function(
    catalogue,
    n_catalogues = 10,
    format = c("sigverse", "matrix"),
    prefix = "catalogue",
    seed = NULL
) {
  # Assertions
  sigshared::assert_catalogue(catalogue)
  format <- rlang::arg_match(format)

  # Derive parameters
  total_mutations <- sum(catalogue[["count"]])
  probs_per_channel <- sigshared::compute_fraction(catalogue[["count"]])

  # Bootstrap with multinomial sampling
  bootstraps <- sigshared::with_seed(
    seed = seed,
    rmultinom(n = n_catalogues, size = total_mutations, prob = probs_per_channel)
  )

  # Add row/column names and type attributes
  rownames(bootstraps) <- catalogue[["channel"]]
  colnames(bootstraps) <- paste0(prefix, "_", seq_len(n_catalogues))
  attr(bootstraps, "type") <- catalogue[["type"]]

  # Return as matrix
  if (format == "matrix") return(bootstraps)

  # Convert to sigverse-style catalogue list
  sigshared::sig_collection_reformat_matrix_to_list(bootstraps)
}

