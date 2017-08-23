#! /usr/bin/env python
#
# Copyright (C) 2015-17 California Institute of Technology

DESCRIPTION = "ztf_sim: Schedule simulator for the Zwicky Transient Facility"
LONG_DESCRIPTION = """\
Schedule simulator for the Zwicky Tranisent Facility
"""

DISTNAME = 'ztf_sim'
MAINTAINER = 'Eric Bellm'
MAINTAINER_EMAIL = 'ecbellm@uw.edu'
URL = 'https://github.com/ZwickyTransientFacility/ztf_sim/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ZwickyTransientFacility/ztf_sim/'
VERSION = '0.0.2.dev'

try:
    from setuptools import setup
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import numpy
    except ImportError:
        install_requires.append('numpy')
    try:
        import scipy
    except ImportError:
        install_requires.append('scipy')
    try:
        import astropy
    except ImportError:
        install_requires.append('astropy')
    try:
        import astroplan
    except ImportError:
        install_requires.append('astroplan')
    try:
        import pandas
    except ImportError:
        install_requires.append('pandas')
    try:
        import sklearn
    except ImportError:
        install_requires.append('sklearn')
    try:
        import sklearn_pandas
    except ImportError:
        install_requires.append('sklearn_pandas')
    try:
        import xgboost
    except ImportError:
        install_requires.append('xgboost')
    try:
        import transitions
    except ImportError:
        install_requires.append('transitions')
    try:
        import gurobipy
    except ImportError:
        install_requires.append('gurobipy')

    return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    setup(name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        install_requires=install_requires,
        packages=['ztf_sim'],
        classifiers=[
                     'Intended Audience :: Science/Research',
                     'Programming Language :: Python :: 2.7',
                     'License :: OSI Approved :: BSD License',
                     'Topic :: Scientific/Engineering :: Visualization',
                     'Operating System :: POSIX',
                     'Operating System :: Unix',
                     'Operating System :: MacOS'],
          )

