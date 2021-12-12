# Causal discovery with potential outcome

This is a code repository for the final project for STAT256 Causal Discovery at UC Berkeley.

## File structure
- `R/`: R files that implement the causal discovery methods.
- `expr/`: R files that handle the numerical simulations and draw the figures
- `data/`: store the results of the numerical simulations.
- `latex/`: tex files for the report.

## Usage
To reproduce the figures in the report, run
```
expr/Rscript do_expr.R pos_signal.R
```
