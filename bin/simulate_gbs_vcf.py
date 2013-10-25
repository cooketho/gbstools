#!/usr/bin/env python
from gbstools import pyms
import subprocess
import sys
import argparse

USAGE = """
simulate_gbs_vcf.py --draws <number of random draws by ms>
                    --samples <number of diploid samples>
                    -t <scaled mutation rate>
"""

DESCRIPTION = """
Draw random segregating sites under a neutral coalescent model with MS (Hudson, 2002).
Use MS results to simulate variation at restriction sites in GBS data.
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('--draws', dest='draws', type=int, help='number of random draws by ms')
parser.add_argument('--samples', dest='n', type=int, help='number of diploid samples')
parser.add_argument('-t', '--theta', dest='t', type=float, help='scaled mutation rate theta')
parser.add_argument('-i', '--input', dest='input', default=None, help='ms results file to use as input')
parser.add_argument('--frag_len', dest='frag_len', type=int, default=500, help='mean digest fragment length (default=500)')
parser.add_argument('--site_len', dest='site_len', type=int, default=6, help='recognition site length (default=6)')
parser.add_argument('--read_len', dest='read_len', type=int, default=101, help='read length (default=101)')
parser.add_argument('--seeds', dest='seeds', type=int, nargs=3, default=[0, 0, 0], help='ms seeds')
parser.add_argument('--variants_only', dest='variants_only', action="store_true", help='output only variant sites')
parser.add_argument('--sampled_only', dest='sampled_only', action="store_true", help='output only sampled sites')
parser.add_argument('--missing', dest='missing', type=float, default=0.25, help='missingness threshold (default=0.25)')
args = parser.parse_args()

try:
    instream = open(args.input, 'r')
except:
    s = args.seeds
    cmd = "ms %i %i -t %f -seeds %i %i %i " % (args.n * 2, args.draws, 
                                               args.t, s[0], s[1], s[2])
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    instream = p.stdout

outstream = sys.stdout

reader = pyms.Reader(instream, fraglen=args.frag_len, sitelen=args.site_len, 
                     readlen=args.read_len, variantsonly=args.variants_only)
n = reader.samples / 2

analysis = ("%s, "  % reader.command,
            "%s, " % reader.seeds,
            "fraglen=%i, " % args.frag_len,
            "readlen=%i, " % args.read_len,
            "sitelen=%i, " % args.site_len,
            "seeds=%s, " % args.seeds,
            "variants_only=%r, " % args.variants_only,
            "sampled_only=%r, " % args.sampled_only,
            "missing=%f" % args.missing)

sitenum = 0
snp = 0
header = "##fileformat=VCFv4.1\n"
header += "##analysis=" + '\t'.join(analysis) + "\n"
header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
header += '\t'.join([str(i) for i in range(n)]) + '\n'
outstream.write(header)
for site in reader:
    sitenum += 1
    if site.ac and site.dcount < 2 * n:    # Ignore monomorphic sites.
        for i in range(len(site.ac)):
            # Ignore sites where ACgbs=0 if --sampled_only argument is used.
            if args.sampled_only and site.ac_sampled[i] == 0:
                continue
            # Ignore sites that fail the missingness threshold.
            elif site.missing[i] > n * args.missing:
                continue
            else:
                snp += 1
                # Get the genotypes sampled by GBS (may differ from true genotypes).
                genotypes = site.genotypes[i]
                # What SNP alleles is the `-` restriction site allele linked to?
                background = site.background[i]
                if background is None:
                    background = '.'
                # Get the ancestral restriction site allele.
                if 0 in site.minus_haplotype:
                    ancestral_rs = '-'
                    if 1 in site.minus_haplotype:
                        ancestral_rs = '+-'
                    else:
                        ancestral_rs = '-'
                else:
                    ancestral_rs = '+'
                    
                    ancestral_rs = '+'
                info = "DCount=%i;" % site.dcount
                info += "AC=%i;" % site.ac[i]
                info += "ACgbs=%i;" % site.ac_sampled[i]
                info += "Missing=%i;" % site.missing[i]
                info += "Background=%s;" % background
                info += "AncestralRS=%s" % ancestral_rs
                fields = ['1', str(sitenum), str(snp)] + ['.'] * 4 + [info, 'GT'] + genotypes
                output = '\t'.join(fields)
                outstream.write(output + '\n')
