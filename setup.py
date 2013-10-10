from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
    CYTHON = True
except ImportError:
    CYTHON = False

requires = []

try:
    import collections
    collections.OrderedDict
except AttributeError:
    requires.append('ordereddict')

cmdclass = {}
ext_modules = []

if CYTHON:
    ext_modules += [Extension("em", ["gbstools/em.pyx"])]
    cmdclass.update({'build_ext':build_ext})
else:
    ext_modules += [Extension("em", ["gbstools/em.c"])]

setup(
    name='GBStools',
    version='0.1.0',
    author='Tom Cooke',
    author_email='cooketho@gmail.com',
    packages=['gbstools'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    license='LICENSE.txt',
    description='Bioinformatics tools for GBS',
    long_description=open('README.txt').read(),
    install_requires=[
        "PyVCF >= 0.6.3",
        "pysam >= 0.7.5",
        "numpy >= 1.3.0",
        "scipy >= 0.7.0",
    ],
)
