# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

This project uses a pyenv virtualenv named `genonets`. Activate it before running any commands:

```bash
pyenv activate genonets
```

## Commands

```bash
# Install dependencies
poetry install

# Run all tests
pytest genonets/test/

# Run a single test file
pytest genonets/test/test_rna.py

# Run a single test by name
pytest genonets/test/test_rna.py::TestRna::test_robustness

# Lint
pylint genonets/

# Format
black genonets/

# Test coverage
coverage run -m pytest genonets/test/ && coverage report

# Build docs
cd docs && make html
```

## Architecture

Genonets is a framework for building and analyzing **genotype networks** — graphs where nodes are biological sequences (RNA/DNA/Protein/Binary) and edges connect sequences differing by a single mutation.

### Main Entry Points

- **`genonets/interface.py`** (`Genonets` class) — the public API. Typical usage:
  ```python
  gn = Genonets(config_params)
  gn.create()    # build networks
  gn.analyze()   # run analyses
  gn.save()      # write output files
  ```
- **`genonets/cmdl_handler.py`** — CLI wrapper around the same `Genonets` class

### Data Flow

```
Input file (genotype-phenotype map)
  → reader.py (InReader)          # parse genotypes + scores per phenotype
  → graph_utils.py (NetworkBuilder) # create igraph networks, add edges between 1-neighbors
  → analysis_handler.py (AnalysisHandler) # run selected analyses
  → writer.py (Writer)            # write GML network files + TSV result files
```

### Key Modules

| Module | Role |
|--------|------|
| `interface.py` | Public API, owns all major objects |
| `graph_utils.py` | `NetworkBuilder` — creates and manipulates `igraph` networks |
| `analysis_handler.py` | Coordinates analyses, supports parallel processing |
| `seq_bit_impl.py` | Alphabet-specific sequence operations (RNA/DNA/Protein/Binary) |
| `config.py` | Configuration parsing and validation |
| `constants.py` | All enums: `AlphabetTypes`, `AnalysisTypes`, `ErrCodes` |
| `errors.py` | Custom exception hierarchy |

### Analysis Modules

Each analysis type has a dedicated `*_functions.py` file with a corresponding analyzer class:

- `peak_functions.py` — fitness peaks/plateaus
- `path_functions.py` — evolutionary paths
- `epistasis_functions.py` — magnitude/sign/reciprocal epistasis
- `robustness_functions.py` — mutational robustness
- `evolvability_functions.py` — evolvability
- `accessibility_functions.py` — phenotype accessibility
- `overlap_functions.py` — phenotypic overlap
- `structure_functions.py` — network structure (clustering, diameter, etc.)
- `landscape_functions.py` — fitness landscape properties

### Sequence Representation

Sequences are internally represented as **bit strings** for efficient 1-neighbor computation. The `seq_bit_lib.py` (abstract) / `seq_bit_impl.py` (concrete) pair handles alphabet-specific encoding. DNA supports reverse complement handling; all types share the same neighbor-finding interface.

### Test Suite

Tests live in `genonets/test/` and validate against ground-truth data in `genonets/test/data/`. Tests are organized by alphabet type (`test_rna.py`, `test_dna.py`, `test_protein.py`, `test_binary.py`) and by analysis type (`test_peaks.py`, `test_paths.py`, `test_epistasis.py`). The package build excludes `genonets/test/` from distribution.

### Sample Scripts

`genonets/sample/` contains runnable examples: `minimal.py`, `simple.py`, `custom.py`, `parallel.py`, `selective.py` — useful references for expected API usage patterns.
