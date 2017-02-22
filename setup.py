#! /usr/bin/env python
#
# Copyright (C) 2015 California Institute of Technology

DESCRIPTION = "ztf_sim: Schedule simulator for the Zwicky Transient Facility"
LONG_DESCRIPTION = """\
Schedule simulator for the Zwicky Tranisent Facility
"""

DISTNAME = 'ztf_sim'
MAINTAINER = 'Eric Bellm'
MAINTAINER_EMAIL = 'ebellm@caltech.edu'
URL = 'https://github.com/ZwickyTransientFacility/ztf_sim/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ZwickyTransientFacility/ztf_sim/'
VERSION = '0.0.1.dev'

try:
    from setuptools import setup
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

if __name__ == "__main__":

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
        install_requires=['numpy', 'scipy', 'matplotlib', 'pandas'],
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

