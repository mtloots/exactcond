"""Build the shared exact-conditional C back-end at INSTALL time, as an Extension.

The earlier setup.py subclassed setuptools' build_py and compiled with distutils' new_compiler.
distutils is gone from Python 3.12 and current setuptools no longer supports that subclassing, so
pip failed with "type object 'build_py' has no attribute '__getattr__'" on every platform (12 Oct
2026). Declaring the sources as an Extension hands the compile to pip, inside build isolation, and
lets real platform wheels be built, after which users need no compiler at all. Same pattern as
arcstat-py.

init.c is excluded deliberately: it is R's symbol-registration glue and references R headers.
"""
import glob
import os

from setuptools import Extension, setup

# Paths must be RELATIVE to setup.py: setuptools rejects absolute source paths.
SRC = os.path.join("src", "exactcond")

sources = sorted(
    s for s in glob.glob(os.path.join(SRC, "*.c"))
    if os.path.basename(s) != "init.c"
)

# Leading underscore: it carries no PyInit and is loaded with ctypes.CDLL by _core.py.
ext = Extension(
    name="exactcond._libexactcond",
    sources=sources,
    include_dirs=[SRC],
    libraries=[] if os.name == "nt" else ["m"],
    optional=False,
)

setup(ext_modules=[ext])
