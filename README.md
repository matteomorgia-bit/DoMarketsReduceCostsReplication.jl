# Replicating "Do Markets Reduce Costs?"

This repository contains a Julia replication project for Fabrizio, Rose, and Wolfram (2007), "Do Markets Reduce Costs? Assessing the Impact of Regulatory Restructuring on US Electric Generation Efficiency."

The project is part of Florian Oswald's Computational Economics course.

## Target Outputs

The replication focuses on:

- Figure 1
- Figure 2
- Table 3
- Table 4
- Table 5

## Links

- Original paper: https://www.aeaweb.org/articles?id=10.1257/aer.97.4.1250
- Replication package: https://www.openicpsr.org/openicpsr/project/116286/version/V1/view
- Course website: https://floswald.github.io/CompEcon/

## Repository Structure

- `src/`: reusable Julia code
- `scripts/`: scripts used to reproduce figures and tables
- `data/raw/`: original raw data files, not modified
- `data/processed/`: cleaned intermediate data
- `output/figures/`: generated figures
- `output/tables/`: generated tables
- `images/`: images included in the report
- `test/`: tests for the Julia code
- `report.qmd`: main replication report

## How to Run

This project uses Julia. From the repository root, activate the project environment with:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

