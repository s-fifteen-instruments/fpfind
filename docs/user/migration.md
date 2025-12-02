# Migration

## v2 to v3

Some options have been deprecated in favor of newer syntax:

1. `-v`/`--verbosity` has been replaced with `-p`/`--quiet`, with a default log level of DEBUG instead of WARNING.
   - `-v2`/`-vv` and above is equivalent to not setting any verbosity flags (DEBUG)
   - `-v1`/`-v` is equivalent to `-p` (INFO)
   - No verbosity flags (or `-v0`) is equivalent to `-pp` (WARNING)

1. `--num-wraps-limit` has been removed.

1. Peak threshold of zero `-S0` now has a special meaning.
   - Peaks will be identified by the proximity of peaks found in early and late segments (relying on low probability of similar peaks from noise), and peak significance will be ignored.

1. `--convergence-rate` has been replaced with `-Q`/`--convergence-order`.
   - Convergence rate was from a misinterpretation of uncertainty scaling between iterations, and behaves badly when deviating from the usual values 4 ~ 6. It also unintuitively scales inversely with the reduction factor.
   - Convergence order is a more [technically-correct term](https://en.wikipedia.org/wiki/Rate_of_convergence), and scales proportionally with the reduction factor without any dependence on acquisition time `Ta` and separation time `Ts`.

Behaviour of fpfind has been slightly modified:

1. Duration overflow is now fixed, by reducing the bin width or number of bins during cross-correlation. This toggle is documented in `fpfind.lib.utils.CoarseHistogramStrategy`.

1. Threshold conditions between iterations now correctly check against previous dt iteration values (when not performing coarse search), and ignore no-change (i.e. cannot be resolved).
