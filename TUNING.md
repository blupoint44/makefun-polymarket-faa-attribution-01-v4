# Tuning notes

All numbers below are mean pinball loss (lower is better) under the
exact scorer, on the six official practice cases unless stated
otherwise.

## Headline result

- House baseline aggregate: **0.363589**
- This entry aggregate: **0.269483**
- Per case: p00 0.356103, p01 0.332848, p02 0.386182,
  p03 0.278749, p04 0.205907, p05 0.113054

## Main settings (partition cases)

- Checkpoints at 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 19, 23, 27, 32,
  38, 45, 53, 62, 73, 85, 99, 115, 133, 154, 178, 205, 240 minutes,
  with log-time interpolation between them.
- History window 4,320 minutes; sample weights decay with a
  720-minute half-life; at most 900 change pairs per checkpoint,
  taken on a regular stride.
- Spread shrinkage 0.65 around the empirical median change; outer
  quantiles (tau 0.05 and 0.95) widened by 1.35x, inner quantiles
  left at 1.0x.
- Regime scale: summed mean absolute one-minute move over 60 minutes
  divided by the same over 1,440 minutes, clipped to [0.7, 1.5].
- Quoted levels snapped to prices traded within the last 1,440
  minutes; exact ties resolve to the higher traded price.
- Partition conditioning: active only with at least 8 markets and a
  last-price sum in [0.6, 1.6]; the gap to `n_markets / 24` is
  returned proportionally to price with a 180-minute phase-in.

## Small-family path

Cases with fewer than 8 markets do not behave like slices of the
24-market partition, and the partition settings overfit them. The
separate small-family settings are:

- width (spread shrinkage) **0.80** instead of 0.65,
- tail multiplier **1.6x** on the outer quantiles instead of 1.35x,
- median weight **0.5**: half of each offset is the empirical median
  change itself, rather than none of it.

On pseudo-cases cut from the FAA tape at 60- and 120-minute
horizons, these small-family settings cut pinball loss by **5–9%**
relative to running the same cases with the partition settings,
which is why the two paths are kept separate instead of averaged
into one compromise set of constants.

## Fallback

Markets with fewer than 60 history points, and any market whose main
forecast raises, use the house-style Gaussian fan: sigma
`0.0025 * sqrt(h)` around the last price (0.5 if there is no price at
all), with the normal z-scores -1.644854, -0.67449, 0, 0.67449,
1.644854 for the five taus. Fallback markets are excluded from the
partition conditioning push.
