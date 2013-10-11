import sys
from optparse import OptionParser
import gbstools

# Parse command line arguments from user.
USAGE = """
score_snps.py -i <input vcf file>
              -o <output vcf file> (default=stdout)
              -d <dispersion> (default=2.5)
              -b <tab-delimited file of sample/bam file pairs> (default=None)
              --normfactors <file output from calc_normfactors.py> (default=None)
              --dpmode <use DP data only?>
              --ped <PLINK-format PED file> (default=None. Use for pedigree mode.)
"""

parser = OptionParser(USAGE)
parser.add_option('-i', '--input', dest='i', help='input bam file')
parser.add_option('-o', '--output', dest='o', default=None, help='output bam file')
parser.add_option('-d', '--dispersion', dest='disp', default=2.5, type='float', help='output bam file')
parser.add_option('-b', '--bam', dest='bamlist', default=None, help='list of sample/bam file pairs')
parser.add_option('-n', '--normfactors', dest='nf', default=None, help='normalization factors file')
parser.add_option('-s', '--samples', dest='samples', default=None, help='samples to include')
parser.add_option('--dpmode',dest='dpmode', action="store_true", help='use DP data only')
parser.add_option('--ped',dest='ped', default=None, help='PED file for nuclear family')
(opt, args) = parser.parse_args()

if opt.o:
   outstream = open(opt.o, 'w')
else:
   outstream = sys.stdout

reader = gbstools.Reader(filename=opt.i, bamlist=opt.bamlist, norm=opt.nf, 
                         disp=opt.disp, ped=opt.ped, samples=opt.samples,
                         dpmode=opt.dpmode)
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
