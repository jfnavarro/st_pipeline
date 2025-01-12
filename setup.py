from setuptools import setup, Extension  # type: ignore
from Cython.Build import cythonize  # type: ignore
import numpy as np

extensions = [
    Extension(
        "stpipeline.common.cdistance",
        ["stpipeline/common/cdistance.pyx"],
        include_dirs=[np.get_include()],
    ),
]

setup(
    name="stpipeline",
    version="2.0.0",
    packages=["stpipeline", "stpipeline.common"],
    ext_modules=cythonize(extensions, language_level="3"),
    include_package_data=True,
    zip_safe=False,
)
