========
GBStools
========

GBStools provides bioinformatics tools for genotyping-by-sequencing data, 
including allele frequency estimation at variant restriction sites.

Installation
============
tar -xvzf GBStools-X.X.X.tar.gz
cd gbstools
python setup.py build
python setup.py install


Quick start
===========
For a quick introduction see docs/TUTORIAL.txt
For the GBStools model, see docs/gbstools_notes.pdf
For more details about using gbstools, see below::


Performing likelihood ratio tests for restriction site variants
===============================================================

To use GBStools to identify restriction site variants that could cause 
genotyping errors, first create a ``Reader`` object that parses a VCF file of 
SNPs you are interested in and then read in the first SNP in the VCF::

       >>> import gbstools
       >>> reader = gbstools.Reader(open('gbstools/test/sim.vcf', 'r'))
       >>> snp = reader.next()

Next, you can calculate maximum likelihood estimates of the frequency of the 
non-cutter restriction site allele by expectation-maximization (EM)::

       >>> while not snp.check_convergence(snp.param['H1']):
       >>>     snp.update_param(snp.param['H1'])
       >>> print snp.param['H1'][-1]
       {'fail': False, 'phi': [0.949, 0.050, 0.00021], 'lambda': 40.008, 'loglik': -39.19, 'delta': 0.0}

Here ``phi`` contains estimates of the frequency of the reference allele, 
alternative allele, and the non-cutter restriction site allele. The ``lambda`` 
variable is an estimate of what the mean coverage would be if all the 
chromosomes were sampled (i.e. no allelic dropout). The ``delta`` variable is 
an estimate of the restriction digest failure rate. See docs/gbstools_notes.pdf 
for more details on the model. The ``fail`` flag is set to True if the boundary 
conditions on the estimates are not met, or something else goes wrong with the 
calculation. And ``loglik`` is the base-e log likelihood of the data given the
current parameters.

You can estimate the null-hypothesis parameters the same way (This time the 
non-cutter allele is restricted to frequency 0).

       >>> while not snp.check_convergence(snp.param['H0']):
       >>>     snp.update_param(snp.param['H0'])
       >>> print snp.param['H0'][-1]
       {'fail': False, 'phi': [0.949, 0.050, 0.0], 'lambda': 40.0, 'loglik': -39.285, 'delta': 0.0}

Now that we know the log-likelihood under H0 and H1, we can calculate the 
likelihood ratio.

       >>> print snp.likelihood_ratio()
       0.178548400541

If you want to create a new vcf that contains the likelihood ratio and other data
in the INFO field, first create a ``Writer`` object. GBStools will copy header 
data from a template vcf (in this case copied from the ``Reader`` object)::

       >>> writer = gbstools.Writer('test.vcf', template=reader)
       >>> snp.lik_ratio = snp.likelihood_ratio()
       >>> snp.update_info()
       >>> writer.write_record(snp)

In most cases your sample libraries will each contain a different number of 
reads, and so you need to provide GBStools with a set of normalization factors
so that the coverage can be compared across samples in a meaningful way in order
to identify which ones are lower than expected (i.e. carry a non-cutter 
restriction site allele). 

       >>> reader = gbstools.Reader(open('gbstools/test/sim.vcf', 'r'), norm='gbstools/test/normfactors.txt')

The format of the normfactors file is:

       insert_size    sample_1    sample_2    ...
       100            1.1         0.9         ...
       101            0.8         1.2         ...
       .	      .		  .	      .
       .	      .		  .	      .

Each row contains the counts of inserts of a given size for each sample, 
normalized to 1.0. See docs/TUTORIAL for instructions on how to make this file 
from your BAM files. For single-end GBS reads the insert data is infered 
from the mapping location and in-silico digest of the reference genome. For 
paired-end GBS data it is extracted from the BAM file. For RAD-seq data, insert
sizes are random, so this file contains only a single row with insert size 
``NA``. If no file is specified, the normalization factors default to 1.0.
In VCF files output by ``Writer``, the normalization factor and insert sizes
are listed in the FORMAT field as ``NF`` and ``INS`` respectively.

GBStools uses information from the SAMPLE fields of the VCF, including the PL
and DP sample fields. If your vcf doesn't have these, or you would like to
calculate these values directly from the alignments, you can pass the ``Reader``
object a file listing a set of bam files::

       >>> reader = gbstools.Reader(open('gbstools/test/sim.vcf', 'r'), bamlist='gbstools/test/bamlist.txt')

The format of the bamlist file is below. The bam files must be indexed.
  
       sample_1    sample_1.bam
       sample_2    sample_2.bam
       .	   .
       .	   .

Running GBStools in pedigree-mode
=================================

By default GBStools assumes genotypes are drawn from an infinite population that
is in Hardy-Weinberg equilibrium. If your SNP data comes from a nuclear familiy 
(i.e. two parents and their offspring), you can change the default model to 
pedigree-mode by passing in a PLINK-formatted PED file::

       >>> reader = gbstools.Reader(open('gbstools/test/sim.ped.vcf', 'r'), ped='gbstools/test/sim.ped')

This EM runs considerably slower than the default EM at the moment, because a 
separate EM is actually run for each possible combination of parental genotypes.
This is a recent addition to GBStools, so after some optimization the runtime 
should improve. If you are interested in this feature, check for updates later.
