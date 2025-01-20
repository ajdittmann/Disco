#!/usr/bin/env python3
import site
from setuptools import setup, Extension
import numpy as np
import imp


# PEP 517 Workaround for edittable user installs
site.ENABLE_USER_SITE = True

version = imp.load_source('discopy.version', 'version.py')


calc_dir = "../Calc"
calcInc = [np.get_include()]
calcInc.append(calc_dir)

calcSources = ["discopy/calc/calcmodule.c",
                calc_dir+"/bondi.c",
                calc_dir+"/integrate.c",
                calc_dir+"/magnetosonic.c"]
calcDeps = [calc_dir+"/calc.h",
            calc_dir+"/integrate.h"]

calcModule = Extension("discopy.calc", sources=calcSources, 
                        include_dirs=calcInc, depends=calcDeps)

setup(
    name='discopy',
    version=version.version,
    description='Python Interface and Analysis Tools for the DISCO code.',
    packages=['discopy'],
    ext_modules=[calcModule],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C",
        "Intended Audience :: Science/Research"],
    install_requires=['numpy>=1.11','h5py>=2.7'],
    extras_require={'docs': ['numpydoc']}
)
