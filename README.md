# fpfind

An add-on implementation of frequency compensation for qcrypto.

Requirements:

* Python 3.8 and above, running in Linux
* `gcc`, if running `freqcd` (preferably in PATH for auto-compilation)

Install the library (alternatively, clone and install locally with `pip install -e.`):

```bash
pip3 install git+https://github.com/s-fifteen-instruments/fpfind.git
```

Binaries and scripts will be exposed to the path; commonly used scripts are listed below.

```bash
fpfind -t {TIMESTAMPS1} -T {TIMESTAMPS2}
freqcd -X -udF fpipe -f 568 < {TIMESTAMPS}
[costream -V5 ... |] freqservo -V5 -udF fpipe -f 568
parse_timestamps -A1 -X -p {TIMESTAMPS}
```

----

## Limitations

* The FFT buffer size is limited to 2**31 bins due to the implicit casting to `int32` performed internally by `np.bincount` on older versions of `numpy`. This corresponds to a buffer order value upper bounded to `q = 31` for `fpfind`. To bypass this limitation, supply an alternative implementation for `np.bincount`.

## Troubleshooting

Certain issues may appear when attempting an install on RaspbianOS:

* Importing `numpy` yields the error message stating `libopenblas.so` could not be found (this is the underlying linear algebra library for Numpy); installing the `libopenblas-dev` library fixes this, e.g. `apt install libopenblas-dev`
