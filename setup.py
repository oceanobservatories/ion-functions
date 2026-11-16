"""
Build configuration for the two optional C extensions. All other project
metadata lives in pyproject.toml. If an extension fails to compile (e.g. no
C compiler), the build continues without it; only the QC functions and
ion_functions.data.polycals depend on these extensions.
"""
import sys

import numpy as np
from setuptools import Extension, setup

include_dirs = [np.get_include(), "extensions/"]
libraries = [] if sys.platform == "win32" else ["m"]

extensions = [
    Extension(
        "ion_functions.qc.qc_extensions",
        sources=["ion_functions/qc/qc_extensions.pyx", "extensions/stuck.c", "extensions/spike.c",
                 "extensions/gradient.c", "extensions/utils.c", "extensions/time_utils.c"],
        include_dirs=include_dirs,
        libraries=libraries,
        optional=True,
    ),
    Extension(
        "ion_functions.data.polycals",
        sources=["ion_functions/data/polycals.pyx", "extensions/polycals.c"],
        include_dirs=include_dirs,
        libraries=libraries,
        optional=True,
    ),
]

setup(ext_modules=extensions)
