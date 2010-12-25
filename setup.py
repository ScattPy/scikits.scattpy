#! /usr/bin/env python
# Last Change: Sun Dec 19 11:00 AM 2010 J

# Copyright (C) 2008 Alexander Vinokurov <alexander.a.vinokurov@gmail.com>

descr   = """ScattPy package.

ScattPy provides numerical methods for solving light scattering problem
by non-spherical particles.
"""

import os
import sys

DISTNAME            = 'scikits.scattpy'
DESCRIPTION         = 'ScattPy package'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'Alexander Vinokurov',
MAINTAINER_EMAIL    = 'alexander.a.vinokurov@gmail.com',
URL                 = 'http://vinokurov.github.com/scattpy'
LICENSE             = 'BSD'
DOWNLOAD_URL        = URL
VERSION             = '0.1-git'

import setuptools
from numpy.distutils.core import setup

def configuration(parent_package='', top_path=None, package_name=DISTNAME):
    if os.path.exists('MANIFEST'): os.remove('MANIFEST')

    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name, parent_package, top_path,
                           version = VERSION,
                           maintainer  = MAINTAINER,
                           maintainer_email = MAINTAINER_EMAIL,
                           description = DESCRIPTION,
                           license = LICENSE,
                           url = URL,
                           download_url = DOWNLOAD_URL,
                           long_description = LONG_DESCRIPTION)

    config.set_options(
        ignore_setup_xxx_py = True,
        assume_default_configuration = True,
        delegate_options_to_subpackages = True,
        quiet = True,
        )

#    config.add_subpackage("scikits")
#    config.add_data_files("scikits/__init__.py")

#    config.add_extension('f_utils',
#                         sources=[os.path.join('src', 'f_utils.for')]
#                         )

    config.add_extension('f_utils',
                         sources=[os.path.join('src', 'f_utils.for')]
                         )
    config.add_extension('f_coords',
                         sources=[os.path.join('src', 'f_coords.for')]
                         )

    return config

if __name__ == "__main__":
    setup(configuration = configuration,
        install_requires = 'numpy',
        namespace_packages = ['scikits'],
        packages = setuptools.find_packages(),
        include_package_data = True,
        #test_suite="tester", # for python setup.py test
        zip_safe = True, # the package can run out of an .egg file
        classifiers =
            [ 'Development Status :: 1 - Planning',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering'])
