# Reproducing the SDY2354 analysis

The `mbmm` package does not redistribute participant-level SDY2354 data.
Obtain the study from ImmPort under accession SDY2354 and comply with the
current ImmPort User Agreement, including its acknowledgment and
non-identification requirements.

Prepare a CSV with exactly these analysis columns:

- `subject`: the ImmPort subject pseudonym;
- `arm`: one of `Placebo`, `LAIV8-IIV5`, `LAIV8-IIV5_Adj`, or `IIV8_Adj`;
- `target`: one of `SinaiIgG`, `H18_IgG`, `H2_IgG`, or `H9_IgG`;
- `day`: study day as an integer;
- `value`: a finite positive serum IgG measurement.

The two original placebo groups are pooled as `Placebo`. `SinaiIgG` denotes
the cH6/1 stalk construct used in the manuscript. The preparation step must
preserve the measurement units used by the manuscript analysis.

Run the published-design analysis with:

```sh
MBMM_SDY2354_CSV=/absolute/path/to/sdy2354_igg_long.csv \
  Rscript sdy2354_table3.R
```

The script fixes the candidate design before scoring, fits a separate outer
bank for every held-out participant, restricts each leave-one-out bank to the
two prespecified class regions, fits the Firth gates, and evaluates both tests
with 4,000 treatment-label permutations.
