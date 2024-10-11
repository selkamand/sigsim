test_that("sig_simulate_catalogues_from_signature works", {
  example_sig <- data.frame(
    channel = c("A[A->G]G", "A[A->G]C", "A[A->G]T", "A[A->G]A"),
    type = rep("A>G", 4),
    fraction = c(0.4, 0.1, 0.3, 0.2)
  )

  # Runs without error
  collection = expect_no_error(sig_simulate_catalogues_from_signature(example_sig, n = 100, n_catalogues = 20, seed = 123))

  # Returns a valid catalogue collection
  expect_no_error(sigshared::assert_catalogue_collection(collection))

  # Length is expected
  expect_length(collection, n = 20)

  # All simulated catalogues have 100 mutations in total
  count_sums = unlist(lapply(collection, \(x){sum(x[["count"]])}))
  expect_true(all(count_sums == 100))

  # Rerunning with same random seed produces identical result
  collection2 = sig_simulate_catalogues_from_signature(example_sig, n = 100, n_catalogues = 20, seed = 123)
  expect_identical(collection, collection2)

  # Rerunning with a different seed produces different results
  collection3 = sig_simulate_catalogues_from_signature(example_sig, n = 100, n_catalogues = 20, seed=1000)
  expect_false(identical(collection2, collection3))

  # format = matrix argument runs without error
  mx <- expect_no_error(sig_simulate_catalogues_from_signature(example_sig, n = 100, n_catalogues = 20, format = "matrix", seed = 123))

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


test_that("sig_simulate_catalogues_from_signatures works", {
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
  collection = expect_no_error(sig_simulate_catalogues_from_signatures(sig_collection, model, n = 100, n_catalogues = 20, seed = 123))

  # Returns a valid catalogue collection
  expect_no_error(sigshared::assert_catalogue_collection(collection))

  # Length is expected
  expect_length(collection, n = 20)

  # All simulated catalogues have 100 mutations in total
  count_sums = unlist(lapply(collection, \(x){sum(x[["count"]])}))
  expect_true(all(count_sums == 100))

  # Rerunning with same random seed produces identical result
  collection2 = sig_simulate_catalogues_from_signatures(sig_collection, model, n = 100, n_catalogues = 20, seed = 123)
  expect_identical(collection, collection2)

  # Rerunning with a different seed produces different results
  collection3 = sig_simulate_catalogues_from_signatures(sig_collection, model, n = 100, n_catalogues = 20, seed = 1000)
  expect_false(identical(collection2, collection3))

  # format = matrix argument runs without error
  mx <- expect_no_error(sig_simulate_catalogues_from_signatures(sig_collection, model, format = "matrix", n = 100, n_catalogues = 20, seed = 123))

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
