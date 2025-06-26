test_that("sig_simulate_single works", {
  example_sig <- data.frame(
    channel = c("A[A->G]G", "A[A->G]C", "A[A->G]T", "A[A->G]A"),
    type = rep("A>G", 4),
    fraction = c(0.4, 0.1, 0.3, 0.2)
  )

  # Runs without error
  collection = expect_no_error(sig_simulate_single(example_sig, n = 100, n_catalogues = 20, seed = 123))

  # Returns a valid catalogue collection
  expect_no_error(sigshared::assert_catalogue_collection(collection))

  # Length is expected
  expect_length(collection, n = 20)

  # All simulated catalogues have 100 mutations in total
  count_sums = unlist(lapply(collection, \(x){sum(x[["count"]])}))
  expect_true(all(count_sums == 100))

  # Rerunning with same random seed produces identical result
  collection2 = sig_simulate_single(example_sig, n = 100, n_catalogues = 20, seed = 123)
  expect_identical(collection, collection2)

  # Rerunning with a different seed produces different results
  collection3 = sig_simulate_single(example_sig, n = 100, n_catalogues = 20, seed=1000)
  expect_false(identical(collection2, collection3))

  # format = matrix argument runs without error
  mx <- expect_no_error(sig_simulate_single(example_sig, n = 100, n_catalogues = 20, format = "matrix", seed = 123))

  # Returns a matrix
  expect_true(is.matrix(mx))

  # Matrix has the expected dimensions (ncol = n_catalogues & nrow = number of channels in sig )
  expect_true(ncol(mx) == 20)
  expect_true(nrow(mx) == nrow(example_sig))

  # Matrix columns sum to 100 mutations
  expect_true(all(colSums(mx) == 100))

  # Rownames are channels (and in the right order)
  expect_identical(rownames(mx), example_sig$channel)
})


test_that("sig_simulate_mixed works", {
  sig_collection <- list(
    sig1 = data.frame(
      channel = c("A[A->G]G", "A[A->G]C", "A[A->G]T"),
      type = rep("A>G", 3L),
      fraction = c(0.4, 0.1, 0.5)
    ),
    sig2 = data.frame(
      channel = c("A[A->G]G", "A[A->G]C", "A[A->G]T"),
      type = rep("A>G", 3L),
      fraction = c(0.4, 0.1, 0.5)
    )
  )

  model <- c(sig1 = 0.3, sig2 = 0.7)

  # Runs without error
  collection = expect_no_error(sig_simulate_mixed(sig_collection, model, n = 100, n_catalogues = 20, seed = 123))

  # Returns a valid catalogue collection
  expect_no_error(sigshared::assert_catalogue_collection(collection))

  # Length is expected
  expect_length(collection, n = 20)

  # All simulated catalogues have 100 mutations in total
  count_sums = unlist(lapply(collection, \(x){sum(x[["count"]])}))
  expect_true(all(count_sums == 100))

  # Rerunning with same random seed produces identical result
  collection2 = sig_simulate_mixed(sig_collection, model, n = 100, n_catalogues = 20, seed = 123)
  expect_identical(collection, collection2)

  # Rerunning with a different seed produces different results
  collection3 = sig_simulate_mixed(sig_collection, model, n = 100, n_catalogues = 20, seed = 1000)
  expect_false(identical(collection2, collection3))

  # format = matrix argument runs without error
  mx <- expect_no_error(sig_simulate_mixed(sig_collection, model, format = "matrix", n = 100, n_catalogues = 20, seed = 123))

  # Returns a matrix
  expect_true(is.matrix(mx))

  # Matrix has the expected dimensions (ncol = n_catalogues & nrow = number of channels in sig )
  expect_true(ncol(mx) == 20)
  expect_true(nrow(mx) == nrow(sig_collection$sig1))

  # Matrix columns sum to 100 mutations
  expect_true(all(colSums(mx) == 100))

  # Rownames are channels (and in the right order)
  expect_true(setequal(rownames(mx), sig_collection$sig1$channel))
})

test_that("sig_simulate_perfect works", {
  sigs <- list(
    sig1 = data.frame(channel = c("C[A>G]C", "C[A>G]T"), type = "A>G", fraction = c(0.7, 0.3)),
    sig2 = data.frame(channel = c("C[A>G]C", "C[A>G]T"), type = "A>G", fraction = c(0.2, 0.8))
  )
  model <- c(sig1 = 0.25, sig2 = 0.75)
  n <- 200

  # Test sigverse format
  perfect <- expect_no_error(sig_simulate_perfect(sigs, model, n = n, format = "sigverse"))

  # Sum of counts = n
  expect_equal(sum(perfect$count), n)

  # Fraction sums to 1
  expect_equal(sum(perfect$fraction), 1, tolerance = 1e-6)

  # Correct columns
  expect_true(all(c("channel", "type", "count", "fraction") %in% names(perfect)))

  # Channels match expected
  expect_setequal(perfect$channel, c("C[A>G]C", "C[A>G]T"))

  # Test matrix format
  mat <- sig_simulate_perfect(sigs, model, n = n, format = "matrix")
  expect_true(is.matrix(mat))
  expect_equal(sum(mat), n)
  expect_setequal(rownames(mat), c("C[A>G]C", "C[A>G]T"))
})

test_that("sig_simulate_stratified works", {
  sigs <- list(
    sig1 = data.frame(channel = c("A[A>G]T", "A[A>G]C"), type = "A>G", fraction = c(0.5, 0.5)),
    sig2 = data.frame(channel = c("A[A>G]T", "A[A>G]C"), type = "A>G", fraction = c(0.6, 0.4))
  )
  model <- c(sig1 = 0.3, sig2 = 0.7)
  n <- 100
  n_cat <- 5

  result <- expect_no_error(sig_simulate_stratified(sigs, model, n = n, n_catalogues = n_cat, seed = 42))

  # Correct structure
  expect_named(result, c("mutations_per_signature", "catalogues"))

  # Check mutation counts per signature
  expect_equal(result$mutations_per_signature$sig1, round(0.3 * n))
  expect_equal(result$mutations_per_signature$sig2, round(0.7 * n))

  # Each catalogue sums to n
  total_counts <- vapply(result$catalogues, function(cat) sum(cat$count), numeric(1))
  expect_true(all(total_counts == n))

  # Is catalogue collection
  expect_no_error(sigshared::assert_catalogue_collection(result$catalogues))

  # Test matrix format
  result_matrix <- sig_simulate_stratified(sigs, model, n = n, n_catalogues = n_cat, format = "matrix", seed = 42)
  expect_true(is.matrix(result_matrix$catalogues))
  expect_equal(ncol(result_matrix$catalogues), n_cat)
  expect_equal(sum(result_matrix$catalogues), n * n_cat)
})


test_that("sig_simulate_single handles edge cases", {
  example_sig <- data.frame(
    channel = c("A[A->G]G", "A[A->G]C"),
    type = rep("A>G", 2),
    fraction = c(0.6, 0.4)
  )

  # n = 0 (no mutations should be sampled)
  out_zero <- sig_simulate_single(example_sig, n = 0, n_catalogues = 2, seed = 1)
  expect_true(all(vapply(out_zero, \(x) sum(x$count), numeric(1)) == 0))

  # n_catalogues <= 0 should throw an error empty list
  expect_error(sig_simulate_single(example_sig, n = 100, n_catalogues = 0), "must be greater than `0`")

  # Probabilities don't sum to 1 (should fail explicitly)
  sig_bad_frac <- data.frame(channel = c("A", "B"), type = "X", fraction = c(10, 30))
  expect_error(sig_simulate_single(sig_bad_frac, n = 100, n_catalogues = 1), "NOT a valid signature")

  # Negative n → should throw error
  expect_error(sig_simulate_single(example_sig, n = -5, n_catalogues = 1), "n must be greater than or equal to `0`")
})


# Test Bootstraps ---------------------------------------------------------

test_that("sig_bootstrap_catalogue returns correct number of catalogues", {
  catalogue <- sigshared::example_catalogue()
  boot <- sig_bootstrap_catalogue(catalogue, n_catalogues = 5)
  expect_type(boot, "list")
  expect_length(boot, 5)
})

test_that("sig_bootstrap_catalogue preserves total mutation count", {
  catalogue <- sigshared::example_catalogue()
  total <- sum(catalogue$count)
  boot <- sig_bootstrap_catalogue(catalogue, n_catalogues = 3)
  for (b in boot) {
    expect_equal(sum(b$count), total)
  }
})

test_that("sig_bootstrap_catalogue returns correct format = matrix", {
  catalogue <- sigshared::example_catalogue()
  mat <- sig_bootstrap_catalogue(catalogue, n_catalogues = 4, format = "matrix")
  expect_true(is.matrix(mat))
  expect_equal(ncol(mat), 4)
  expect_equal(nrow(mat), nrow(catalogue))
  expect_equal(sum(mat[, 1]), sum(catalogue$count))
})

test_that("sig_bootstrap_catalogue sets reproducible seed", {
  catalogue <- sigshared::example_catalogue()
  boot1 <- sig_bootstrap_catalogue(catalogue, n_catalogues = 3, seed = 42)
  boot2 <- sig_bootstrap_catalogue(catalogue, n_catalogues = 3, seed = 42)
  expect_equal(boot1, boot2)
})

test_that("sig_bootstrap_catalogue names output correctly", {
  catalogue <- sigshared::example_catalogue()
  boot <- sig_bootstrap_catalogue(catalogue, n_catalogues = 2, prefix = "test")
  expect_named(boot, c("test_1", "test_2"))
})

