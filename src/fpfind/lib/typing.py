import os
from typing import Iterator, List, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias

from fpfind import NP_PRECISEFLOAT

Number: TypeAlias = Union[int, float, np.integer, np.floating]
Integer: TypeAlias = Union[int, np.integer]
Float: TypeAlias = Union[float, np.floating]
Complex_: TypeAlias = Union[np.complex64, np.complex128, np.complex256]
Complex: TypeAlias = Union[complex, Complex_]

NumberArray: TypeAlias = Union[NDArray[np.number], List[Number]]
IntegerArray: TypeAlias = Union[NDArray[np.integer], List[Integer]]
FloatArray: TypeAlias = Union[NDArray[np.floating], List[Float]]
ComplexArray: TypeAlias = Union[NDArray[Complex_], List[Complex]]

TimestampArray: TypeAlias = NDArray[Union[NP_PRECISEFLOAT, np.uint64]]
DetectorArray: TypeAlias = NDArray[np.uint32]
TimestampDetectorArray: TypeAlias = Tuple[TimestampArray, DetectorArray]
TimestampDetectorStream: TypeAlias = Iterator[TimestampDetectorArray]

PathLike: TypeAlias = Union[os.PathLike, str]
Mask: TypeAlias = NDArray[np.bool_]
