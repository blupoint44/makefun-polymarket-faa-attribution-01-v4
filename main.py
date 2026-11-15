"""Quantile forecaster for a family of related probability markets.

The idea is deliberately plain: recent one-minute moves, reweighted so
that fresh observations matter more than stale ones, give an empirical
distribution of price changes at a handful of look-ahead distances.
Interpolating those distributions in log-time yields a full horizon,
the whole picture is widened or narrowed by a short-versus-long
volatility ratio, and finally every quoted level is pulled back onto
a price the market has actually traded at recently.  Markets that
belong to a partition of one underlying event also get a gentle push
so that their last prices, which should sum to a known total, stop
drifting away from it.

Only the standard library is used and nothing here is random, so a
given input always produces the same output.
"""

import bisect
import glob
import json
import math
import os

LEVELS_TAU = (0.05, 0.25, 0.5, 0.75, 0.95)
NORMAL_Z = (-1.644854, -0.67449, 0.0, 0.67449, 1.644854)
FALLBACK_SIGMA = 0.0025

# Look-ahead distances, in minutes, at which the empirical change
# distribution is evaluated directly.  Everything between them is
# interpolated, so they are spaced roughly geometrically.
CHECKPOINTS = (
    1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 19, 23, 27, 32, 38, 45, 53,
    62, 73, 85, 99, 115, 133, 154, 178, 205, 240,
)

HISTORY_CAP = 4320        # minutes of history fed to the change sampler
DECAY_HALF_LIFE = 720.0   # minutes for a sample weight to halve
SPREAD_FACTOR = 0.65      # shrinkage applied to empirical spreads
TAIL_FACTORS = (1.35, 1.0, 1.0, 1.0, 1.35)
PAIR_BUDGET = 900         # most change samples used per checkpoint
SNAP_LOOKBACK = 1440      # window whose traded prices form the grid
VOL_SHORT = 60
VOL_LONG = 1440
VOL_RATIO_LO = 0.7
VOL_RATIO_HI = 1.5
FAMILY_SIZE = 24          # markets in a complete partition
REVERSION_TAU = 180.0     # minutes for the partition push to phase in
PARTITION_MIN = 8         # fewest markets that still look like a family

# A case with very few markets is not a slice of the big partition, so
# it is treated more cautiously: wider bands, fatter tails, and half
# of the forecast anchored on the empirical median change.
SMALL_SPREAD = 0.80
SMALL_TAILS = (1.6, 1.0, 1.0, 1.0, 1.6)
SMALL_MEDIAN_WEIGHT = 0.5


def tidy_price(value):
    """Clamp a probability into range and round it to six decimals."""
    bounded = min(1.0, max(0.0, float(value)))
    return round(bounded * 1000000.0) / 1000000.0 + 0.0


def fallback_rows(centre, horizon):
    """Gaussian fan around the last price, used when history is thin."""
    rows = []
    for step in range(1, horizon + 1):
        spread = FALLBACK_SIGMA * math.sqrt(step)
        rows.append([tidy_price(centre + z * spread) for z in NORMAL_Z])
    return rows


def interpolated_quantiles(ordered_values, ordered_weights):
    """Read the five target quantiles off a weighted empirical sample.

    The sample arrives sorted by value with a weight attached to each
    entry.  Walking the cumulative weight up to each tau fraction of
    the total and interpolating linearly inside the straddling entry
    gives the corresponding quantile.
    """
    count = len(ordered_values)
    total = math.fsum(ordered_weights)
    if count == 0 or total <= 0.0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    answers = []
    for tau in LEVELS_TAU:
        goal = tau * total
        pos = 0
        seen = 0.0
        while pos < count - 1 and seen + ordered_weights[pos] < goal:
            seen += ordered_weights[pos]
            pos += 1
        weight_here = ordered_weights[pos]
        if weight_here <= 0.0:
            share = 0.0
        else:
            share = (goal - seen) / weight_here
            share = min(1.0, max(0.0, share))
        following = ordered_values[pos + 1] if pos + 1 < count else ordered_values[pos]
        answers.append(ordered_values[pos] + share * (following - ordered_values[pos]))
    return answers


def checkpoint_offsets(series, regime, spread=SPREAD_FACTOR,
                       tails=TAIL_FACTORS, median_weight=0.0):
    """Empirical change quantiles at every checkpoint distance.

    Pairs are formed between prices a checkpoint apart.  Each pair is
    weighted by how recent its ending price is, with the weight
    halving every half-life, and when there are more pairs than the
    budget allows a regular stride keeps the sample spread over the
    whole history instead of just its tail.
    """
    size = len(series)
    shrink = math.exp(math.log(0.5) / DECAY_HALF_LIFE)
    age_weights = [shrink ** exponent for exponent in range(size)]
    table = []
    for distance in CHECKPOINTS:
        starts = size - distance
        if starts <= 0:
            table.append([0.0, 0.0, 0.0, 0.0, 0.0])
            continue
        stride = max(1, starts // PAIR_BUDGET)
        # Stepping by the stride visits exactly the start indices
        # 0, stride, 2*stride, ... below starts, in order.
        sampled = []
        for index in range(0, starts, stride):
            change = series[index + distance] - series[index]
            weight = age_weights[size - 1 - (index + distance)]
            sampled.append((change, weight))
        sampled.sort(key=lambda pair: pair[0])
        quantiles = interpolated_quantiles(
            [pair[0] for pair in sampled], [pair[1] for pair in sampled])
        middle = quantiles[2]
        row = []
        for slot in range(5):
            centred = (quantiles[slot] - middle) * spread * tails[slot] * regime
            row.append(middle * median_weight + centred)
        table.append(row)
    return table


def closest_traded(grid, target):
    """Return the traded price nearest to target, preferring the
    higher neighbour whenever target sits exactly in between."""
    if not grid:
        return target
    place = bisect.bisect_left(grid, target)
    if place <= 0:
        return grid[0]
    if place >= len(grid):
        return grid[-1]
    below = grid[place - 1]
    above = grid[place]
    if above - target <= target - below:
        return above
    return below


def mean_abs_step(series, window):
    """Average absolute one-minute move over the recent window."""
    usable = min(window, len(series) - 1)
    if usable <= 0:
        return 0.0
    total = 0.0
    for idx in range(len(series) - usable, len(series)):
        total += abs(series[idx] - series[idx - 1])
    return total / usable


def current_regime(all_series):
    """Short-horizon activity relative to the long-horizon norm."""
    recent_total = 0.0
    usual_total = 0.0
    for series in all_series:
        recent_total += mean_abs_step(series, VOL_SHORT)
        usual_total += mean_abs_step(series, VOL_LONG)
    if usual_total <= 0.0 or recent_total <= 0.0:
        return 1.0
    ratio = recent_total / usual_total
    return min(VOL_RATIO_HI, max(VOL_RATIO_LO, ratio))


def stretch_offsets(table, horizon):
    """Expand checkpoint offsets to every minute of the horizon.

    Distances before the first checkpoint reuse it and distances past
    the last one reuse that; in between, blending is done against the
    logarithm of the distance because the checkpoints themselves are
    spaced geometrically.
    """
    result = []
    log_points = [math.log(point) for point in CHECKPOINTS]
    for step in range(1, horizon + 1):
        if step <= CHECKPOINTS[0]:
            result.append(list(table[0]))
            continue
        if step >= CHECKPOINTS[-1]:
            result.append(list(table[-1]))
            continue
        seg = 0
        while CHECKPOINTS[seg + 1] < step:
            seg += 1
        fraction = ((math.log(step) - log_points[seg])
                    / (log_points[seg + 1] - log_points[seg]))
        blended = []
        for slot in range(5):
            low = table[seg][slot]
            blended.append(low + fraction * (table[seg + 1][slot] - low))
        result.append(blended)
    return result


def enforce_order(row):
    """Make a quantile row non-decreasing in place and return it."""
    for slot in range(1, len(row)):
        if row[slot] < row[slot - 1]:
            row[slot] = row[slot - 1]
    return row


def single_market(series, horizon, regime, small_family):
    last = series[-1]
    if len(series) < 60:
        return fallback_rows(last, horizon)
    grid = sorted(set(series[-SNAP_LOOKBACK:]))
    recent = series[-HISTORY_CAP:]
    if small_family:
        offsets = checkpoint_offsets(recent, regime, spread=SMALL_SPREAD,
                                     tails=SMALL_TAILS,
                                     median_weight=SMALL_MEDIAN_WEIGHT)
    else:
        offsets = checkpoint_offsets(recent, regime)
    stretched = stretch_offsets(offsets, horizon)
    finished = []
    for offsets_row in stretched:
        quoted = []
        for offset in offsets_row:
            snapped = closest_traded(grid, tidy_price(last + offset))
            quoted.append(tidy_price(snapped))
        finished.append(enforce_order(quoted))
    return finished


def nudge_partition(all_rows, all_series, is_fallback, horizon):
    """Shift a partition family so its implied total drifts home.

    If the last prices of the whole family add up to more (or less)
    than a complete partition should, each trustworthy market is
    moved by a share of the discrepancy proportional to its own
    price, phased in over a few hours and snapped back onto traded
    prices so the adjustment stays on the market's own grid.
    """
    lasts = [series[-1] for series in all_series]
    total_now = math.fsum(lasts)
    discrepancy = total_now - len(all_series) / FAMILY_SIZE
    if discrepancy == 0.0:
        return
    trusted_total = 0.0
    for flag, value in zip(is_fallback, lasts):
        if not flag:
            trusted_total += value
    if trusted_total <= 0.0:
        return
    for position in range(len(all_series)):
        if is_fallback[position]:
            continue
        proportion = lasts[position] / trusted_total
        if proportion == 0.0:
            continue
        grid = sorted(set(all_series[position][-SNAP_LOOKBACK:]))
        rows = all_rows[position]
        for step in range(1, horizon + 1):
            shift = -(1.0 - math.exp(-step / REVERSION_TAU)) * discrepancy * proportion
            if shift == 0.0:
                continue
            moved = [tidy_price(closest_traded(grid, tidy_price(v + shift)))
                     for v in rows[step - 1]]
            rows[step - 1] = enforce_order(moved)


def handle_case(payload):
    horizon = int(payload.get("horizon_minutes", 240))
    markets = payload.get("markets", [])
    histories = []
    for market in markets:
        histories.append([float(v) for v in market.get("prices", [])])
    try:
        regime = current_regime(histories)
    except Exception:
        regime = 1.0
    small_family = len(markets) < PARTITION_MIN
    forecasts = []
    fallback_flags = []
    for market, history in zip(markets, histories):
        used_fallback = len(history) < 60
        try:
            rows = single_market(history, horizon, regime, small_family)
        except Exception:
            centre = history[-1] if history else 0.5
            rows = fallback_rows(centre, horizon)
            used_fallback = True
        fallback_flags.append(used_fallback)
        forecasts.append({"id": market.get("id"), "q": rows})
    # A market with no history at all contributes nothing to the
    # family total; every other market contributes its last price.
    last_sum = math.fsum(h[-1] if h else 0.0 for h in histories)
    if len(markets) >= PARTITION_MIN and 0.6 <= last_sum <= 1.6:
        try:
            plain_rows = [entry["q"] for entry in forecasts]
            nudge_partition(plain_rows, histories, fallback_flags, horizon)
        except Exception:
            pass
    return {"forecasts": forecasts}


def main():
    source_dir = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
    target_dir = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")
    os.makedirs(target_dir, exist_ok=True)
    names = sorted(glob.glob(os.path.join(source_dir, "*.in")))
    for path in names:
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        outcome = handle_case(payload)
        stem = os.path.basename(path)
        if stem.endswith(".in"):
            stem = stem[:-3]
        destination = os.path.join(target_dir, stem + ".out")
        with open(destination, "w", encoding="utf-8") as handle:
            json.dump(outcome, handle)


if __name__ == "__main__":
    main()
