import pyms
import subprocess
from optparse import OptionParser

"""Draw random segregating sites under a neutral coalescent model with MS (Hudson, 2002).
Store MS output in a file and use these results to simulate variation at restriction
sites in GBS data.
"""

USAGE = """
simulate_geno.py -N <number of sites>
                 -n <number of samples>
                 -t <scaled mutation rate>
                 -p <output prefix>
                 --frag_len <mean digest fragment length> (default=500)
                 --site_len <enzyme recognition site length> (default=6)
                 --read_len <read length> (default=101)
                 --seeds <3 comma-separated ms random seeds> (default=0,0,0)
                 --variants_only <output only restriction site variants?>
"""

parser = OptionParser(USAGE)
parser.add_option('-N', '--sites', dest='N', type="int", help='number of sites')
parser.add_option('-n', '--samples', dest='n', type="int", help='number of samples')
parser.add_option('-t', '--theta', dest='t', type="float", help='scaled mutation rate theta')
parser.add_option('-p', '--prefix', dest='prefix', type="string", help='output prefix')
parser.add_option('--frag_len', dest='frag_len', type="int", default=500, help='mean digest fragment length')
parser.add_option('--site_len', dest='site_len', type="int", default=6, help='recognition site length')
parser.add_option('--read_len', dest='read_len', type="int", default=101, help='read length')
parser.add_option('--seeds', dest='seeds', type="string", default="0,0,0", help='ms seeds')
parser.add_option('--variants_only', dest='variants_only', action="store_true", help='output only variant sites')
(opt, args) = parser.parse_args()

s = opt.seeds.split(',')
cmd = "ms %i %i -t %f -seeds %s %s %s " % (opt.n * 2, opt.N, opt.t, s[0], s[1], s[2])
cmd = cmd + "> %s.ms" % (opt.prefix)
subprocess.call(cmd, shell=True)

reader = pyms.Reader(open(opt.prefix + '.ms', 'r'), N=opt.frag_len, L=opt.site_len, 
                     readlen=opt.read_len, variantsonly=opt.variants_only)

sitenum = 0
snp = 0
nmissing_thresh = 0.2
outfile = open(opt.prefix + '.vcf', 'w')
header = "##fileformat=VCFv4.1\n"
header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
header += '\t'.join([str(i) for i in range(opt.n)])
outfile.write(header + '\n')
for site in reader:
    sitenum += 1
    if site.ac and site.dcount < 2 * opt.n:    # Ignore monomorphic sites.
        for i in range(len(site.ac)):
            if site.ac_sampled[i] > 0 and site.missing[i] <= opt.n * nmissing_thresh:
                snp += 1
                genotypes = site.genotypes[i]
                background = site.background[i]
                if background is None:
                    background = '.'
                info = "DCount=%i;" % site.dcount
                info += "AC=%i;" % site.ac[i]
                info += "ACgbs=%i;" % site.ac_sampled[i]
                info += "Missing=%i;" % site.missing[i]
                info += "Background=%s" % background
                fields = ['1', str(sitenum), str(snp)] + ['.'] * 4 + [info, 'GT'] + genotypes
                output = '\t'.join(fields)
                outfile.write(output + '\n')
