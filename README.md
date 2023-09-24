
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigsim <img src="man/figures/sigsim_hex_96dpi.png" align="right" height="120" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigsim)](https://CRAN.R-project.org/package=sigsim)
<!-- badges: end -->

**sigsim** simplifies simulation of mutational signature decompositons.

**sigsim** is part of the
[**sigverse**](https://github.com/selkamand/sigverse)

sigsim is in early development and not yet ready for use.

## Philosophy

Simulating datasets from mutational signature collections is essential
for benchmarking the performance of mutational signature analysis tools,
particularly identifying failure modes of tool.

For example, In cancer contexts we’re often interested in a single base
substitution known as ‘Signature 3’. This is because the signature seems
to correlate with an underlying biological phenomena that relates a
patient’s likelihood of responding to certain therapies.

For important signatures like this, its important we understand how well
we can actually identify the signature, and what other signatures may
cause false positives / false negatives. These sorts of questions can be
addressed using simulation experiments.

This package allows you to create simulated datasets, where by ‘dataset’
we mean per-sample decompositions describing a tally of ‘mutations’
matching each ‘channel’ in a signature.

To generate these datasets we first paramaterise a multinomial
distribution based on a single mutational signature OR some combination
of multiple mutational signatures. We then randomly sample this
multimonium distribution ‘n’ times (where n = number of mutations we
want to simulate for a sample).

This can be done with: `sig_simulate_decomposition()`

When you simulate decompositions by random sampling a multinomial model
in this way, the randomness of the sampling adds noise. For some
experiments, you may want to include decompositions that represent
perfect (‘noiseless’) combinations of two signatures.

You can create these, using the `noise = FALSE` argument of
`sig_simulate_decomposition()`

## Installation

You can install the development version of sigsim from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/sigsim")
```

## Quick Start

``` r
library(sigsim) # Signature simulation
library(sigstash) # Signature collections
library(sigvis) # Visualisation of signatures


# Load signature collections from sigstash
signatures = sig_load('COSMIC_V3.3.1_SBS')


# Randomly simulate a 400 mutations, by random sampling SBS1
sig_simulate_decomposition(signatures, model = c('SBS1' = 1), n = 400)

# Randomly simulate a 400 mutations, where ~30% come from SBS1 and ~70% come from SBS3
sig_simulate_decomposition(signatures, model = c('SBS1' = 0.3, 'SBS3' = 0.7), n = 400)

# Simulate a decomposition with 400 mutations, which represents a perfect addition of 30% SBS1 and 70% SBS3
sig_simulate_decomposition(signatures, model = c('SBS1' = 0.3, 'SBS3' = 0.7), n = 400, noise = FALSE)
```