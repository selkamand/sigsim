
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigsim <img src="man/figures/sigsim_hex_96dpi.png" align="right" height="120" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigsim)](https://CRAN.R-project.org/package=sigsim)
<!-- badges: end -->

**sigsim** simplifies simulation of mutational signature catalogues.

**sigsim** is part of the
[**sigverse**](https://github.com/selkamand/sigverse)

sigsim is in early development and not yet ready for use.

## Philosophy

Simulating datasets from mutational signature collections is essential
for benchmarking the performance of mutational signature analysis tools,
particularly identifying failure modes.

For example, In cancer contexts we’re often interested in a single base
substitution known as ‘Signature 3’. This is because the signature seems
to correlate with an underlying biological phenomena that relates a
patient’s likelihood of responding to certain therapies.

For important signatures like this, its important to understand how well
we can actually identify the signature, and what other signatures may
cause false positives / false negatives. These questions can be
addressed using simulation experiments.

**sigsim** simplifies creation of simulated datasets, where by ‘dataset’
we mean per-sample catalogues describing a tally of ‘mutations’ matching
each ‘channel’ in a signature.

To generate these datasets we first paramaterise a multinomial
distribution based on a single mutational signature OR some combination
of multiple mutational signatures. We then randomly sample this
multimonium distribution ‘n’ times (where n = number of mutations we
want to simulate for a sample).

When you simulate catalogues by random sampling a multinomial model in
this way, the randomness of the sampling adds noise. For some
experiments, you may want to include catalogues that represent perfect
(‘noiseless’) combinations of two signatures.

You can create these ‘perfect’ datasets using the `noise = FALSE`
argument of `sig_simulate_catalogue()`

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
sig_simulate_catalogue(signatures, model = c('SBS1' = 1), n = 400)

# Randomly simulate a 400 mutations, where ~30% come from SBS1 and ~70% come from SBS3
sig_simulate_catalogue(signatures, model = c('SBS1' = 0.3, 'SBS3' = 0.7), n = 400)

# 'Simulate' a catalogue with 400 mutations, which represents a perfect addition of 30% SBS1 and 70% SBS3
sig_simulate_catalogue(signatures, model = c('SBS1' = 0.3, 'SBS3' = 0.7), n = 400, noise = FALSE)
```

### Simulating Datasets for Sensitivity, Inter-Signature Interference, and Overfitting Assessment

The datasets simulated above are well-suited for evaluating
**sensitivity** and **inter-signature interference**. However, if your
focus is on **assessing overfitting**, you may need more control over
the level and type of ‘noise’ added to the model.

In the previous examples, noise is introduced via random sampling from
the multinomial distribution. But for overfitting assessment, we might
want this noise to diverge more significantly from the signatures we’re
modeling.

#### Adding Controlled Noise to Simulations

One approach is to simulate **n1** mutations that follow a ‘perfect’
signature and then add **n2** mutations from a clearly defined ‘noise’
model. This raises the question: **what should our noise model be?**

Here are two options:

1.  **Uniform Distribution**  
    A simple choice is to generate noise by randomly distributing the
    **n2** mutations across all available channels. This represents a
    completely uniform noise model, where mutations are allocated
    without reference to any specific signature.

2.  **Low-Cosine-Similarity Signature**  
    A more sophisticated approach is to simulate noise using an
    artificially created signature with a low cosine similarity to
    **all** known signatures in a given collection. This simulates the
    addition of noise from a mutational process not represented in the
    current dataset.  
    This approach is particularly useful when applying mutational
    signatures in novel contexts—such as using signature databases
    developed for adult cancers in childhood cancer research.
    Differences in mutational processes between disease types may lead
    to signatures that don’t fully capture the true spectrum of
    mutational signatures, increasing the risk of overfitting.

#### Simulating Overfitting with Dropped Signatures

Instead of artificially creating a low-similarity signature, another
effective strategy is to simulate mutation catalogues using a
combination of known signatures, but then **drop one signature** when
performing the fitting. If the fitting algorithm attempts to explain the
missing signature by over-relying on other available signatures, it
indicates overfitting.

This method allows you to test how well the algorithm can handle
‘unknown’ components—i.e., real mutational processes that aren’t
represented in your signature collection.
