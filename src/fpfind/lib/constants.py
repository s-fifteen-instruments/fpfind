import enum
import warnings

import numpy as np


class TSRES(enum.Enum):
    """Stores timestamp resolution information.

    Values assigned correspond to the number of units within a
    span of 1 nanosecond.
    """

    NS2 = 0.5  # S-Fifteen TDC1 timestamp
    NS1 = 1
    PS125 = 8  # CQT Red timestamp
    PS4 = 256  # S-Fifteen TDC2 timestamp


EPOCH_LENGTH = 1 << 29  # from filespec
FCORR_AMAXBITS = -13  # from 'freqcd.c'
NTP_MAXDELAY_NS = 200e6  # for very *very* asymmetric channels

# Derived constants
EPOCH_DURATION = EPOCH_LENGTH * 1e-9  # in seconds
MAX_FCORR = 2**FCORR_AMAXBITS

# Compilations of numpy that do not include support for 128-bit floats will not
# expose 'np.float128'. We map such instances directly into a 64-bit float instead.
# Note that some variants implicitly map 'np.float128' to 'np.float64' as well.
#
# The availability of this alias depends on the available build flags for numpy,
# and is platform-dependent. This is usually extended-precision (80-bits) padded
# to 128-bits for efficiency.
NP_PRECISEFLOAT = np.float64
if hasattr(np, "float128"):
    NP_PRECISEFLOAT = np.float128
else:
    warnings.warn(
        "Extended-precision floats unsupported in current numpy build on this "
        "platform: falling back to 64-bit floats instead.\n\n"
        "Hint: To avoid precision loss when using time taggers with <125ps "
        "precision, read timestamps as 64-bit integers instead, e.g. "
        "'read_a1(..., fractional=False)'."
    )


class PeakFindingFailed(ValueError):
    def __init__(
        self,
        message,
        significance=None,
        resolution=None,
        dt1=None,
        dt2=None,
        dt=None,
        df=None,
    ):
        self.message = message
        self.s = significance
        self.r = resolution
        self.dt1 = dt1
        self.dt2 = dt2
        self.dt = dt
        self.df = df
        super().__init__(message)

    def __str__(self):
        text = self.message
        suppl = []
        if self.s is not None:
            suppl.append(f"S={self.s:8.3f}")
        if self.r is not None:
            suppl.append(f"r={self.r:5.0f}")
        if self.dt1 is not None:
            suppl.append(f"dt1={self.dt1:11.0f}")
        if self.dt2 is not None:
            suppl.append(f"dt2={self.dt2:11.0f}")
        if self.dt is not None:
            suppl.append(f"dt={self.dt:11.0f}")
        if self.df is not None:
            suppl.append(f"df={self.df * 1e6:.4f}ppm")
        if suppl:
            text = f"{text} ({', '.join(suppl)})"
        return text


class PeakFindingTerminated(RuntimeError):
    pass
