#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''CNV analysis tools'''

setup(
    name='cnv_analysis',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['cnv_analysis'],
    package_dir={'cnv_analysis': 'src'},
    entry_points={
        'console_scripts': ['collapse_cnvs = cnv_analysis.collapse_cnvs:main']
    },
    url='https://github.com/bjpop/cnv_analysis',
    license='LICENSE',
    description=('CNV analysis tools'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["networkx"],
)