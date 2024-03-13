# Pairwise Alignment Algorithms

This repository contains all my implementations (in Rust) for the algorithms for pairwise alignment covered in the course *Computational Challenges in Bioinformatics* at Ghent University.

## Installation

Install Rust using the information provided [here](https://www.rust-lang.org/tools/install).

## Compilation

Run following command in the root of the repository to compile all the executables.

```sh
cargo build --release
```

## Overview

This repository contains the code for the following algorithms:

- [Needleman-Wunsch](./needleman_wunsch)
- [Hirschberg](./hirschberg)
- [Banded global alignment](./banded_global_alignment)
- [Needleman-Wunsch using delta encoded scoring](./needleman_wunsch_delta_encoded)
- [BitPAl](./bitpal)

### Executing BitPAl
The actual implementation of the BitPAl algorithm can be found in [`bitpal/src/lib.rs`](bitpal/src/lib.rs).

Computing the pairwise global alignment score using the BitPAl algorithm can be done in following steps (after building the executables):
1) Navigate to the release directory ([`target/release`](target/release)) starting at the root of this repository.
2) Execute following command:
```shell
./bitpal -f <input_file.fasta>
```

