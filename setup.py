from distutils.core import setup
from distutils.extension import Extension

requires = []

try:
    import collections
    collections.OrderedDict
except AttributeError:
    requires.append('ordereddict')

try:
    import vcf
except ImportError:
    requires.append('PyVCF')

try:
    import pysam
except ImportError:
    requires.append('pysam')

try:
    import numpy
except ImportError:
    requires.append('numpy')

try:
    import scipy
except ImportError:
    requires.append('scipy')

try:
    from Cython.Distutils import build_ext
    CYTHON = True
except ImportError:
    CYTHON = False

cmdclass = {}
ext_modules = []

if CYTHON:
    ext_modules += [Extension('em', ['gbstools/em.pyx'])]
    cmdclass.update({'build_ext':build_ext})
else:
    ext_modules += [Extension('em', ['gbstools/em.c'])]

scripts = ['bin/annotate_pe_bam.py',
           'bin/annotate_se_bam.py',
           'bin/digest_to_bed.py',
           'bin/likelihood_ratio.py',
           'bin/make_gbstools_bed.py',
           'bin/mapping_summary.py',
           'bin/normfactors.py',
           'bin/simulate_dp.py',
           'bin/simulate_fasta.py',
           'bin/simulate_fastq.py',
           'bin/simulate_gbs_vcf.py',
           'bin/simulate_ped_vcf.py']

setup(
    name='GBStools',
    version='0.1.0',
    author='Tom Cooke',
    author_email='cooketho@gmail.com',
    packages=['gbstools'],
    scripts=scripts,
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    license='LICENSE.txt',
    description='Bioinformatics tools for GBS',
    long_description=open('README.txt').read(),
    requires=requires,
)
