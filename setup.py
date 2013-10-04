from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("em", ["gbstools/em.pyx"])]

setup(
    name='GBStools',
    version='0.1.0',
    author='Tom Cooke',
    author_email='cooketho@gmail.com',
    packages=['gbstools'],
    cmdclass={'build_ext':build_ext},
    ext_modules=ext_modules,
    license='LICENSE.txt',
    description='Bioinformatics tools for GBS',
    long_description=open('README.txt').read(),
    install_requires=[
        "PyVCF >= 0.6.3",
        "pysam >= 0.7.5",
        "numpy >= 1.3.0",
        "scipy >= 0.7.0",
        "Cython >= 0.17.4"
    ],
)
