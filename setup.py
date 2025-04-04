#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    print("setuptools not found, falling back to distutils")
    from distutils.core import setup

setup(
    name='pycorn',
    version='0.20',
    author='R. Jaepel',
    packages=['pycorn'],
    requires=["xmltodict"],
    extras_require={'plotting':  ["matplotlib"], 'xlsx-output': ['xlsxwriter'], "processing": ["numpy", "pandas"],
                    "testing": ["pytest"]},
    scripts=['examplescripts/pycorn-bin.py'],
    platforms=['Linux', 'Windows', 'MacOSX'],
    zip_safe=False,
    classifiers=["License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                 "Environment :: Console",
                 "Intended Audience :: Science/Research",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 3.9",],
    package_data={'pycorn': ['docs/*.*']},
    license='GNU General Public License v2 (GPLv2)',
    description='A script to extract data from UNICORN result (res) files',
    long_description=open('README.md').read(),
    url='https://github.com/ronald-jaepel/PyCORN',
)
