#!/usr/bin/env python
from gbstools import pyms
import subprocess
import sys
import argparse

USAGE = """
simulate_gbs_vcf.py -N <number of sites>
                    -n <number of samples>
                    -t <scaled mutation rate>
"""

DESCRIPTION = """
Draw random segregating sites under a neutral coalescent model with MS (Hudson, 2002).
Use MS results to simulate variation at restriction sites in GBS data.
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-N', '--sites', dest='N', type=int, help='number of sites')
parser.add_argument('-n', '--samples', dest='n', type=int, help='number of samples')
parser.add_argument('-t', '--theta', dest='t', type=float, help='scaled mutation rate theta')
parser.add_argument('-p', '--prefix', dest='prefix', default=None, help='output prefix')
parser.add_argument('--frag_len', dest='frag_len', type=int, default=500, help='mean digest fragment length (default=500)')
parser.add_argument('--site_len', dest='site_len', type=int, default=6, help='recognition site length (default=6)')
parser.add_argument('--read_len', dest='read_len', type=int, default=101, help='read length (default=101)')
parser.add_argument('--seeds', dest='seeds', type=int, nargs=3, default=[0, 0, 0], help='ms seeds')
parser.add_argument('--variants_only', dest='variants_only', action="store_true", help='output only variant sites')
parser.add_argument('--missing', dest='missing', type=float, default=0.2, help='missingness threshold (default=0.2)')
args = parser.parse_args()

s = args.seeds
cmd = "ms %i %i -t %f -seeds %i %i %i " % (args.n * 2, args.N, args.t, s[0], s[1], s[2])

try:
    outstream = open(args.prefix + '.vcf', 'w')
    cmd = cmd + "> %s.ms" % (args.prefix)
    subprocess.check_call(cmd, shell=True)
    instream = open(args.prefix + '.ms', 'r')
except:
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    instream = p.stdout
    outstream = sys.stdout

reader = pyms.Reader(instream, N=args.frag_len, L=args.site_len, 
                     readlen=args.read_len, variantsonly=args.variants_only)

sitenum = 0
snp = 0
header = "##fileformat=VCFv4.1\n"
header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
header += '\t'.join([str(i) for i in range(args.n)])
outstream.write(header + '\n')
for site in reader:
    sitenum += 1
    if site.ac and site.dcount < 2 * args.n:    # Ignore monomorphic sites.
        for i in range(len(site.ac)):
            if site.ac_sampled[i] > 0 and site.missing[i] <= args.n * args.missing:
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
                outstream.write(output + '\n')
