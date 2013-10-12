#!/usr/bin/env python
import sys
import argparse
import gbstools

# Parse command line arguments from user.
USAGE = """
polymorphism_test.py -i <input vcf file>
                     -o <output vcf file>
                     -b <tab-delimited file of sample/bam file pairs>
                     --normfactors <output from normfactors.py>
"""

DESCRIPTION = """
Calculate likelihood ratio for absence/presence of restriction site polymorphism
at GBS SNPs with GBStools, and output maximum likelihood estimates of allele
frequency at restriction sites.
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-i', '--input', dest='i', help='input VCF file containing GBS SNPs (sites with >2 alleles will be ignored)', required=True)
parser.add_argument('-o', '--output', dest='o', default=None, help='output VCF file (default is stdout)')
parser.add_argument('-d', '--dispersion', dest='disp', default=2.5, type=float, help='dispersion index used in GBStools EM (default=2.5)')
parser.add_argument('-b', '--bam', dest='bamlist', default=None, help='list of sample/bam file pairs (use bam files instead of VCF to get alignment data)')
parser.add_argument('-n', '--normfactors', dest='nf', default=None, help='normalization factors file used in GBStools EM (produced by normfactors.py)')
parser.add_argument('-s', '--samples', dest='samples', default=None, help='samples to include (excluded samples will remain in output VCF, but will not be used in EM)')
parser.add_argument('--dpmode',dest='dpmode', action="store_true", help='use DP data only; ignore PL data from VCF')
parser.add_argument('--ped',dest='ped', default=None, help='PED file for nuclear family (when specified, pedigree-mode is used)')
args = parser.parse_args()

if args.o:
   outstream = open(args.o, 'w')
else:
   outstream = sys.stdout

reader = gbstools.Reader(filename=args.i, bamlist=args.bamlist, norm=args.nf, 
                         disp=args.disp, ped=args.ped, samples=args.samples,
                         dpmode=args.dpmode)
writer = gbstools.Writer(outstream, template=reader)

for snp in reader:
   if not reader.family:
      while not snp.check_convergence(snp.param['H1']):
         snp.update_param(snp.param['H1'])
         if len(snp.param['H1']) > 40:
            snp.param['H1'][-1]['fail'] = True
      while not snp.check_convergence(snp.param['H0']):
         snp.update_param(snp.param['H0'])
         if len(snp.param['H0']) > 40:
            snp.param['H0'][-1]['fail'] = True
   else:
      for gt in snp.param:
         while not snp.check_convergence(snp.param[gt]):
            snp.update_param(param=snp.param[gt], parental_gt=gt)
            if len(snp.param[gt]) > 40:
               snp.param[gt][-1]['fail'] = True
   snp.lik_ratio = snp.likelihood_ratio()
   snp.update_info()
   writer.write_record(snp)
