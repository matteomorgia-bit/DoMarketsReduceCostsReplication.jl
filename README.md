# Replicating "Do Markets Reduce Costs?"

This repository contains a Julia replication project for Fabrizio, Rose, and Wolfram (2007), "Do Markets Reduce Costs? Assessing the Impact of Regulatory Restructuring on US Electric Generation Efficiency."

The project is part of Florian Oswald's Computational Economics course.

## Replication Targets

The project reproduces:

- Figure 1
- Figure 2
- Tables 3, 4, and 5

The implementation translates the relevant Stata routines from the original replication package into Julia.

## Links

- Original paper: https://www.aeaweb.org/articles?id=10.1257/aer.97.4.1250
- Replication package: https://www.openicpsr.org/openicpsr/project/116286/version/V1/view
- Course website: https://floswald.github.io/CompEcon/
- Online report: https://matteomorgia-bit.github.io/DoMarketsReduceCostsReplication.jl/

## Data

The raw data are not committed to this repository. Download the openICPSR replication package and place the extracted folder here:

```text
data/raw/openicpsr/116286-V1/
```

The scripts expect files such as:

```text
frw1extract_enf.dta
frw1extract_f.dta
inputregs.do
praisiv2.do
```

## Install

From the repository root, start Julia and run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Run Tests

```bash
julia --project=. test/runtests.jl
```

## Run Replication

The package has a single entry point:

```julia
using DoMarketsReduceCostsReplication
DoMarketsReduceCostsReplication.run_all()
```

This runs the data checks, table scripts, and figure scripts.

Individual scripts are also in `scripts/`.

## Report

The main project report is:

```text
report.qmd
```

The rendered online report is available at:

https://matteomorgia-bit.github.io/DoMarketsReduceCostsReplication.jl/

The report is designed to be rendered and deployed with GitHub Pages using the workflow in `.github/workflows/publish.yml`.

