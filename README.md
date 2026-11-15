# polymarket-faa-attribution-01 — quantile forecaster

A deterministic, standard-library-only forecaster for a family of
related probability markets. For every market in a case it emits 240
rows (one per horizon minute, or as many as the case's
`horizon_minutes` asks for) of five quantiles at tau 0.05, 0.25, 0.5,
0.75 and 0.95.

## How it works

1. **Empirical changes, recency-weighted.** At 28 checkpoint distances
   (1 to 240 minutes, roughly geometrically spaced) the program pairs
   prices that far apart over the last 4,320 minutes of history. Each
   pair's weight halves every 720 minutes of age, and at most 900
   pairs per checkpoint are kept, spread by a regular stride.
2. **Shrunk spreads, fattened tails.** The weighted quantiles of those
   changes are centred on their median, shrunk by a factor of 0.65,
   and the two outer quantiles are widened again by 1.35x.
3. **Regime scaling.** The whole fan is scaled by the ratio of the
   mean absolute one-minute move over the last 60 minutes to the same
   measure over 1,440 minutes, summed across the case's markets and
   clipped to [0.7, 1.5].
4. **Log-time interpolation.** Between checkpoints, offsets are
   blended linearly in the logarithm of the horizon minute.
5. **Snapping to traded prices.** Every quoted level is moved to the
   nearest price actually traded in the last 1,440 minutes (ties go to
   the higher price), rounded to six decimals, and each row is forced
   to be non-decreasing.
6. **Partition conditioning.** When a case has at least 8 markets and
   their last prices sum to between 0.6 and 1.6, the family is treated
   as a partition of 24 markets: any gap between the observed sum and
   `n_markets / 24` is pushed back into the non-fallback markets in
   proportion to their prices, phased in with a 180-minute time
   constant and snapped onto each market's traded grid.
7. **Fallbacks.** A market with fewer than 60 history points gets a
   plain Gaussian fan around its last price with sigma
   `0.0025 * sqrt(h)`. Any failure in the main path falls back the
   same way rather than failing the case.

Small cases (fewer than 8 markets) are not treated as partitions:
they use a wider spread (0.80), fatter tails (1.6x on the outer
quantiles) and anchor half of each offset on the empirical median
change. See `TUNING.md` for why.

## Running

```
MAKEFUN_INPUT_DIR=/path/to/input MAKEFUN_OUTPUT_DIR=/path/to/output python3 main.py
```

Every `*.in` file in the input directory is read in sorted order and
a matching `<name>.out` file is written, containing
`{"forecasts": [{"id": ..., "q": [[...] x horizon]}]}`.

No network access, no third-party packages, no randomness: the same
input always yields byte-identical output.

## Validation on the practice kit

Scored with the exact pinball scorer on the six practice cases:

| case | pinball loss |
|------|--------------|
| p00  | 0.356103 |
| p01  | 0.332848 |
| p02  | 0.386182 |
| p03  | 0.278749 |
| p04  | 0.205907 |
| p05  | 0.113054 |
| **aggregate** | **0.269483** |

The house baseline on the same cases aggregates to 0.363589, so this
entry cuts the baseline pinball loss by roughly a quarter.
