"""Build the shared C back-end at install time, so no compiler is needed at import."""
import os, glob, sysconfig
from setuptools import setup
from setuptools.command.build_py import build_py
from distutils.ccompiler import new_compiler

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src", "exactcond")

class Build(build_py):
    def run(self):
        srcs = sorted(f for f in glob.glob(os.path.join(SRC, "*.c"))
                      if os.path.basename(f) != "init.c")
        if srcs:
            cc = new_compiler()
            objs = cc.compile(srcs, output_dir=self.build_temp,
                              extra_preargs=["-O2", "-fPIC"])
            cc.link_shared_lib(objs, "_libexactcond", output_dir=SRC, libraries=["m"])
        super().run()

setup(cmdclass={"build_py": Build}, package_data={"exactcond": ["*.c", "*.so", "*.dylib"]})
