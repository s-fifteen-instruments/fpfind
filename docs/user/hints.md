# Usage hints

## Getting help

Usage help for `fpfind` is provided in increasing verbosity levels:

```bash
fpfind -h    # for basic usage
fpfind -hh   # for fine-tuning compensation
fpfind -hhh  # for testing, changing limits, changing return behaviour
```

## fpfind

Initial test with timestamps or epochs:

```bash
fpfind -t {ALICE} -T {BOB}    # defined per timestamp filespec
fpfind -d {T2DIR} -D {T1DIR}  # defined per qcrypto filespec
```

Optimization of peak finding often requires tuning the initial resolution `-R` (units of ns) and FFT buffer order `-q` (corresponding to the number of FFT bins $N=2^q$). The rate at which the timing uncertainty is reduced between iterations can be adjusted with `-Q` (default of 2.71 corresponds to ~1.44 bits of accuracy per iteration).

```bash
fpfind ... -q23 -R64 -Q3
```

The integrated channel selector can be used to select specific detector channels, to filter out irrelevant channels. Example below uses only channel 1 (`0b0001`) events from reference, and channel 3 and 4 events (`0b1100`) from target.

```bash
fpfind ... -m 1 -M 12
```

Do full frequency pre-compensation scan:

```bash
fpfind ... -P --precomp-step 100e-9 --precomp-ordered --precomp-fullscan
```

Under low signal conditions, the peak significance threshold can be adjusted:

```bash
fpfind ... -S5.6
```

Use liberal peak searching (ignoring peak significance), by setting peak threshold to zero. In practice, useful to reduce the early-late separation as well to minimize peak drift, especially under high frequency offset conditions (e.g. >10ppm).

```bash
fpfind ... -S0 -s1
```

The output of `fpfind` can be customized with `-V`, e.g. `-V5` to only report the compensated frequency, or `-V10` to use `freqcd`-compatible units (2^-34 for frequency, 0.125ps for time).

The correctness of peak finding can be graphically visualized with `fpplot` (currently requires additional dependencies [`kochen`](https://pypi.org/project/kochen/) and [`S15lib`](https://github.com/s-fifteen-instruments/pyS15#) that need to be side-loaded). Swap the program from `fpfind` to `fpplot` directly, and supply the `--df` and `--dt` results from `fpfind` into `fpplot` (assuming the compensating values were calculated, i.e. `-V0`), e.g.

```bash
> fpfind ... -V0
-9.98468151003351e-06   -439017
> fpplot ... -V0 --df=-9.98468151003351e-06 --dt=-439017
```

If `fpfind` needs to be quickly terminated, use the `--interruptible` flag and terminate with `SIGINT`:

```bash
fpfind ... -I
```