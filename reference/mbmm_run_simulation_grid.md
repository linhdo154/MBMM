# Run a manuscript-style simulation grid

Runs independent simulation scenarios sequentially or across multiple
worker processes.

## Usage

``` r
mbmm_run_simulation_grid(grid, productive_bank, n_cores = 1L)
```

## Arguments

- grid:

  A simulation grid returned by \[mbmm_demo_scenario_grid()\].

- productive_bank:

  Predicted productive reference bank.

- n_cores:

  Number of worker processes. Use \`1\` for sequential execution.

## Value

A list containing one fitted result for each row of \`grid\`.
