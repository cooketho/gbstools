#!/usr/bin/env python
import re
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
parser.add_argument('-b', '--bam', dest='bamlist', default=None, help='list of sample/bam file pairs (use bam files instead of VCF to get alignment data)')
parser.add_argument('-n', '--normfactors', dest='nf', default=None, help='normalization factors file used in GBStools EM (default NF=1.0 for all samples). See also normfactors.py.')
parser.add_argument('-s', '--samples', dest='samples', default=None, help='samples to include (excluded samples will remain in output VCF, but will not be used in EM)')
parser.add_argument('--dispersion_slope', dest='disp_slope', default=0.0, type=float, help='slope for linear function of dispersion index vs mean coverage (default=0.0)')
parser.add_argument('--dispersion_intercept', dest='disp_intercept', default=2.5, type=float, help='intecept for linear function of dispersion index vs mean coverage (default=2.5)')
parser.add_argument('--dpmode',dest='dpmode', action="store_true", help='use DP data only; ignore PL data from VCF')
parser.add_argument('--ped',dest='ped', default=None, help='PED file for nuclear family (when specified, pedigree-mode is used)')
parser.add_argument('--intervals',dest='intervals', default=None, type=str, help='samtools-style intervals (e.g. chr1:1-1000)')
parser.add_argument('--debug',dest='debug', action='store_true', help='use debug mode')
args = parser.parse_args()

if args.o:
   outstream = open(args.o, 'w')
else:
   outstream = sys.stdout

reader = gbstools.Reader(filename=args.i, bamlist=args.bamlist, norm=args.nf, 
                         disp_slope=args.disp_slope, disp_intercept=args.disp_intercept,
                         ped=args.ped, samples=args.samples, dpmode=args.dpmode)
                         
writer = gbstools.Writer(outstream, template=reader)

try:
   intervals = re.split('[:-]', args.intervals)
   chrom = intervals.pop(0)
   try:
      start = int(intervals.pop(0))
   except:
      start = None
   try:
      end = int(intervals.pop(0))
   except:
      end = None
   snps = reader.fetch(chrom, start, end)
except:
   snps = (snp for snp in reader)

for snp in snps:
   if not reader.family:
      while not snp.check_convergence(snp.param['H1']):
         param = snp.update_param(snp.param['H1'])
         snp.param['H1'].append(param)
         if len(snp.param['H1']) > 30:
            snp.param['H1'][-1]['fail'] = True
      param = snp.update_param(snp.param['H1'])
      snp.param['H1'][-1]['loglik'] = param['loglik']

      while not snp.check_convergence(snp.param['H0']):
         param = snp.update_param(snp.param['H0'])
         snp.param['H0'].append(param)
         if len(snp.param['H0']) > 30:
            snp.param['H0'][-1]['fail'] = True
      param = snp.update_param(snp.param['H0'])
      snp.param['H0'][-1]['loglik'] = param['loglik']

   else:
      for gt in snp.param:
         while not snp.check_convergence(snp.param[gt]):
            param = snp.update_param(param=snp.param[gt], parental_gt=gt)
            snp.param[gt].append(param)
            if len(snp.param[gt]) > 30:
               snp.param[gt][-1]['fail'] = True
         param = snp.update_param(snp.param[gt])
         snp.param[gt][-1]['loglik'] = param['loglik']

   snp.lik_ratio = snp.likelihood_ratio()
   snp.update_info()

   if args.debug:
      print 'H0:'
      for param in snp.param['H0']:
         print param
      print 'H1:'
      for param in snp.param['H1']:
         print param
      print 'DP: %s' % str([call.DP for call in snp.calls])
      print 'PL: %s' % str([call.PL for call in snp.calls])
      print 'NF: %s' % str([call.NF for call in snp.calls])

   writer.write_record(snp)


