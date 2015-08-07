#!/usr/bin/env python
from gbstools import pyms
import subprocess
import sys
import argparse
from collections import namedtuple

USAGE = """
simulate_gbs_vcf.py --draws <number of random draws by ms>
                    --samples <number of diploid samples>
                    -t <scaled mutation rate>
"""

DESCRIPTION = """
Draw random segregating sites under a neutral coalescent model with MS (Hudson, 2002).
Use MS results to simulate variation at restriction sites in GBS data.
"""

"""INFO fields to be added to the vcf header."""
_Info = namedtuple('Info', ['id', 'num', 'type', 'desc'])
INFO = (_Info('DCount', None, 'Integer', 'Non-cut restriction site allele count'),
        _Info('AC', None, 'Integer', 'Allele count in true genotypes'),
        _Info('ACgbs', None, 'Integer', 'Allele count in genotypes after accounting for GBS allelic dropout'),
        _Info('Hets', None, 'Integer', 'Number of true heterozygous genotypes'),
        _Info('HetsGBS', None, 'Integer', 'Number of apparent heterozygous genotypes after accounting for GBS allelic dropout'),
        _Info('HetsMissing', None, 'Integer', 'Number of true heterozygous genotypes that appear as missing (./.) after accounting for GBS allelic dropout'),
        _Info('Missing', None, 'Integer', 'Number of missing (./.) genotype calls after accounting for GBS allelic dropout'),
        _Info('Background', None, 'String', 'SNP allele(s) which co-occur on same haplotype as non-cut restriction site allele. Possible values: ancestral, derived, both'),
        _Info('AncestralRS', None, 'String', 'Ancestral state of variable restriction site(s). Possible values: + (cut), - (non-cut), +- (multiple variable restriction site alleles with different ancestral states)'))

"""FORMAT fields to be added to the vcf header."""
_Format = namedtuple('Format', ['id', 'num', 'type', 'desc'])
FORMAT = _Format('GT', None, 'String', 'Genotype. Alleles that are unobservable by GBS are denoted by `.` where the true value of the unobservable allele can be either 0 or 1')

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('--draws', dest='draws', type=int, help='number of random draws by ms')
parser.add_argument('--samples', dest='n', type=int, help='number of diploid samples')
parser.add_argument('-t', '--theta', dest='t', type=float, help='scaled mutation rate theta')
parser.add_argument('-i', '--input', dest='input', default=None, help='ms results file to use as input')
parser.add_argument('--frag_len', dest='frag_len', type=int, default=500, help='mean digest fragment length (default=500)')
parser.add_argument('--site_len', dest='site_len', type=int, default=6, help='recognition site length (default=6)')
parser.add_argument('--site_prob', dest='site_prob', type=float, default=0.0011, help='per-bp probability of recognition site')
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
                     siteprob=args.site_prob, readlen=args.read_len, 
                     variantsonly=args.variants_only)
n = reader.samples / 2

analysis = ("%s, "  % reader.command,
            "%s, " % reader.seeds,
            "fraglen=%i, " % args.frag_len,
            "readlen=%i, " % args.read_len,
            "sitelen=%i, " % args.site_len,
            "siteprob=%f, " % args.site_prob,
            "seeds=%s, " % args.seeds,
            "variants_only=%r, " % args.variants_only,
            "sampled_only=%r, " % args.sampled_only,
            "missing=%f" % args.missing)

outstream.write("##fileformat=VCFv4.1\n")
outstream.write("##analysis=" + '\t'.join(analysis) + "\n")
for info in INFO:
    outstream.write("##INFO=<ID=%s,Number=.,Type=%s,Description=\"%s\">\n" % (info.id, info.type, info.desc))
outstream.write("##FORMAT=<ID=%s,Number=.,Type=%s,Description=\"%s\">\n" % (FORMAT.id, FORMAT.type, FORMAT.desc))
header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
header += '\t'.join([str(i) for i in range(n)]) + '\n'
outstream.write(header)
for site in reader:
    if site.ac and site.dcount < 2 * n:    # Ignore monomorphic sites.
        for i in range(len(site.ac)):
            # Ignore sites where ACgbs=0 if --sampled_only argument is used.
            if args.sampled_only and site.ac_sampled[i] == 0:
                continue
            # Ignore sites that fail the missingness threshold.
            elif site.missing[i] > n * args.missing:
                continue
            else:
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
                info += "Hets=%i;" % site.hets[i]
                info += "HetsGBS=%i;" % site.hets_sampled[i]
                info += "Missing=%i;" % site.missing[i]
                info += "Background=%s;" % background
                info += "AncestralRS=%s" % ancestral_rs

                fields = [str(site.sitenum), str(i)] + ['.'] * 5 + [info, 'GT'] + genotypes
                output = '\t'.join(fields)
                outstream.write(output + '\n')
