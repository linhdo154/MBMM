# Simulate one transient or abortive detection pattern

The four patterns reflect the reduced manuscript simulation:

- `noise_like`: no intended signal above the LOD;

- `single_blip`: one intended detectable visit;

- `two_point`: two intended near-LOD detections;

- `fast_rise_drop`: one higher isolated peak with neighboring visits
  remaining near background.

## Usage

``` r
mbmm_simulate_abortive_pattern(
  times,
  lod = 100,
  sigma_obs = 0.5,
  subtype = NULL,
  subtype_prob = c(noise_like = 0.2, single_blip = 0.35, two_point = 0.3, fast_rise_drop
    = 0.15),
  max_log10_above_lod = 1.5,
  target_detect_prob = 0.9,
  seed = NULL
)
```

## Arguments

- times:

  Observation times.

- lod:

  Natural-scale assay limit of detection.

- sigma_obs:

  Observation-noise SD on the log10 scale.

- subtype:

  Optional fixed subtype. When `NULL`, one is sampled.

- subtype_prob:

  Named subtype probabilities.

- max_log10_above_lod:

  Maximum signal above log10(LOD).

- target_detect_prob:

  Probability used to protect intended detections from being lost after
  observation noise.

- seed:

  Optional random seed.

## Value

A data frame containing natural-scale observations, censoring, latent
log10 signal, and subtype.
