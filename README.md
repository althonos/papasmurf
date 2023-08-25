# 🧙‍♂️ PAPASMURF [![Star me](https://img.shields.io/github/stars/althonos/papasmurf.svg?style=social&label=Star&maxAge=3600)](https://github.com/althonos/papasmurf/stargazers)

*A [Platform-Accelerated](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data) Package for Alignment-free [SMURF](https://github.com/NoamShental/SMURF) analysis.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/papasmurf/rust.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/papasmurf/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/papasmurf?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/papasmurf/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![Crate](https://img.shields.io/crates/v/papasmurf.svg?maxAge=600&style=flat-square)](https://crates.io/crates/papasmurf)
[![Docs](https://img.shields.io/docsrs/papasmurf?maxAge=600&style=flat-square)](https://docs.rs/papasmurf)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/papasmurf/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/papasmurf/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/papasmurf.svg?style=flat-square&maxAge=600)](https://github.com/althonos/papasmurf/issues)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/papasmurf/blob/master/CHANGELOG.md)

## 🗺️ Overview

SMURF (Short MUltiple Region Framework) is a method proposed by Fuks 
*et al.*[\[1\]](#ref1) in 2018 for taxonomic profiling of 16S sequencing
data. It uses several PCR-amplified regions inside the 16S rRNA gene to
reach high taxonomic resolution despite the use of short read sequencing.

PAPASMURF is a Rust reimplementation of the SMURF method from scratch. It 
does ***not*** aim at being a 1-to-1 reimplementation of the original 
[MATLAB implementation](https://github.com/NoamShental/SMURF), but allows 
more control over the parameters used in the original to support sequencing 
data of lesser quality.

### 📋 Features

- **primer profiles**: The primers defining each regions are converted into 
  a Position-Specific Scoring Matrix (PSSM) to allow for the fast extraction
  of 16S gene regions when building the database.
- **sparse matrices**: The mapping and reconstruction algorithms are entirely
  implemented with sparse matrices, reducing the memory consumption and 
  accelerating dot-product computations[\[2\]](#ref2).
- **fast k-mer matching**: The full-scan k-mer matching in the mapping phase 
  is implemented using [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data)
  to compare each read to all the k-mers in a database. While of $O(kn)$ 
  runtime complexity, this is actually faster than using a dedicated data
  structure (such as a [trie](https://en.wikipedia.org/wiki/Trie)) to 
  recover k-mers with mismatches.


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/papasmurf/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

<!-- ### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/althonos/papasmurf/blob/master/CONTRIBUTING.md) for more details. -->

## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/papasmurf/blob/master/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.

## ⚖️ License

This library is provided under the open-source
[GPLv3 license](https://choosealicense.com/licenses/gpl-3.0/).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original SMURF authors](https://github.com/NoamShental). It was developed 
by [Martin Larralde](https://github.com/althonos/) during his PhD project at 
the [European Molecular Biology Laboratory](https://www.embl.de/) in the 
[Zeller team](https://github.com/zellerlab).*

## 📚 References

- <a id="ref1">\[1\]</a> Fuks, Garold, Michael Elgart, Amnon Amir, Amit Zeisel, Peter J. Turnbaugh, Yoav Soen, and Noam Shental. ‘Combining 16S RRNA Gene Variable Regions Enables High-Resolution Microbial Community Profiling’. Microbiome 6 (26 January 2018): 17. [doi:10.1186/s40168-017-0396-x](https://doi.org/10.1186/s40168-017-0396-x).
- <a id="ref2">\[2\]</a> Gustavson, Fred G. ‘Two Fast Algorithms for Sparse Matrices: Multiplication and Permuted Transposition’. ACM Transactions on Mathematical Software 4, no. 3 (September 1978): 250–69. [doi:10.1145/355791.355796](https://doi.org/10.1145/355791.355796).