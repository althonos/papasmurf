# 🧙‍♂️ PAPASMURF [![Star me](https://img.shields.io/github/stars/althonos/papasmurf.svg?style=social&label=Star&maxAge=3600)](https://github.com/althonos/papasmurf/stargazers)

*A [Platform-Accelerated](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data) Package for Alignment-free [SMURF](https://github.com/NoamShental/SMURF) analysis.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/papasmurf/python.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/papasmurf/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/papasmurf?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/papasmurf/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Docs](https://img.shields.io/readthedocs/papasmurf/latest?style=flat-square&maxAge=600)](https://papasmurf.readthedocs.io)
[![Crate](https://img.shields.io/crates/v/papasmurf-py.svg?maxAge=600&style=flat-square)](https://crates.io/crates/papasmurf-py)
[![PyPI](https://img.shields.io/pypi/v/papasmurf.svg?style=flat-square&maxAge=600)](https://pypi.org/project/papasmurf)
[![Wheel](https://img.shields.io/pypi/wheel/papasmurf.svg?style=flat-square&maxAge=2678400)](https://pypi.org/project/papasmurf/#files)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/papasmurf?style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/papasmurf)
[![Python Versions](https://img.shields.io/pypi/pyversions/papasmurf.svg?style=flat-square&maxAge=600)](https://pypi.org/project/papasmurf/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/papasmurf.svg?style=flat-square&maxAge=600)](https://pypi.org/project/papasmurf/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/papasmurf/tree/main/papasmurf-py)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/papasmurf/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/papasmurf.svg?style=flat-square&maxAge=600)](https://github.com/althonos/papasmurf/issues)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/papasmurf/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/papasmurf?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/papasmurf)

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

*This is the Python version, there is a [Rust crate](https://crates.io/crates/papasmurf) available as well.*

## 🔧 Installing

<!-- `papasmurf` can be installed directly from [PyPI](https://pypi.org/project/papasmurf/),
which hosts some pre-built wheels for most mainstream platforms, as well as the 
code required to compile from source with Rust:
```console
$ pip install papasmurf
``` -->
<!-- Otherwise, `papasmurf` is also available as a [Bioconda](https://anaconda.org/bioconda/papasmurf)
package:
```console
$ conda install -c bioconda papasmurf
``` -->

In the event you have to compile the package from source, all the required
Rust libraries are vendored in the source distribution, and a Rust compiler
will be setup automatically if there is none on the host machine.


## 💡 Example

Use [Biopython](https://biopython.org) to generate a database from a
file containing 16S gene sequences in FASTA format, for instance the 
[Greengenes database](https://greengenes.secondgenome.com/):

```python
import papasmurf

# Create a database builder with the two given primers
builder = papasmurf.Builder([
    ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"),  # V3-V4 primers
    ("GTGYCAGCMGCCGCGGTAA", "CCGYCAATTYMTTTRAGTTT"), # V4-V5 primers
])

# Extract k-mers from the reference sequences
with gzip.open("gg_13_5.fasta.gz", "rt") as reader:
    for record in Bio.SeqIO.parse(reader, "fasta"):
        builder.add(record.id, str(record.seq))

# Build and index the database
database = builder.to_database()

# Save the database in JSON format
database.dump("gg.json", format="json")
```

Then use the database to map reads from a sample:

```python
# Load database and create a new mapper
database = papasmurf.Database.load("gg.json", format="json")
mapper = papasmurf.Mapper(database)

# Map reads to the k-mers database
with gzip.open("data/Example_L001_R1_001.fastq.gz", "rt") as f1:
    with gzip.open("data/Example_L001_R2_001.fastq.gz", "rt") as f2:
        for r1, r2 in zip(Bio.SeqIO.parse(f1, "fastq"), Bio.SeqIO.parse(f2, "fastq")):
            mapper.add(str(r1.seq), str(r2.seq))
```

Once all the reads have been mapped, compute the final bacterium frequencies:

```python
# Obtain partial mapping result
result = mapper.finish()

# Run the iterative procedure 10 times to estimate the read proportion vector
result.refine(10)

# Print the names of the reference sequences with >5% relative abundance
for (j, name) in enumerate(database.names):
    if result.frequencies[j] > 0.05:
        print(name, result.frequencies[j])
```


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
[Zeller team](https://github.com/zellerlab) with support and testing
from [Fabian Springer](https://github.com/fabispringer).*

*All brand names and product names used in this material are trademarks or registered trademarks of their respective owners. The author/owner is not affiliated with, endorsed by, or sponsored by any product, organization, or company mentioned. Smurf is a registered trademark of Studio Peyo S.A.*

## 📚 References

- <a id="ref1">\[1\]</a> Fuks, Garold, Michael Elgart, Amnon Amir, Amit Zeisel, Peter J. Turnbaugh, Yoav Soen, and Noam Shental. ‘Combining 16S RRNA Gene Variable Regions Enables High-Resolution Microbial Community Profiling’. Microbiome 6 (26 January 2018): 17. [doi:10.1186/s40168-017-0396-x](https://doi.org/10.1186/s40168-017-0396-x).
- <a id="ref2">\[2\]</a> Gustavson, Fred G. ‘Two Fast Algorithms for Sparse Matrices: Multiplication and Permuted Transposition’. ACM Transactions on Mathematical Software 4, no. 3 (September 1978): 250–69. [doi:10.1145/355791.355796](https://doi.org/10.1145/355791.355796).