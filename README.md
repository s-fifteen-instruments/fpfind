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
freqcd -x -df 568 < {TIMESTAMPS}
parse_timestamps -A1 -X -p {TIMESTAMPS}
```
