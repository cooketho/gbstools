#!/usr/bin/env python
import random
import sys
from scipy.stats import nbinom
import math
import numpy
import argparse

USAGE = """
cat <myvcf> | simulate_dp.py --lambda <coverage> > <mynewvcf>
"""

DESCRIPTION = """
Use GT data in a GBS vcf file to simulate AD and PL data, where AD values are
drawn from a negative binomial distribution with mean=lambda and index of
dispersion=d. For calculating PL, the case call error rate is assumed to be
0.001 by default.
"""

parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
parser.add_argument('-l', '--lambda', dest='lamb', type=float, help='mean depth of coverage', required=True)
parser.add_argument('-e', '--epsilon', dest='epsilon', type=float, default=0.001, help='base call error rate (default=0.001)')
parser.add_argument('--dispersion_mean', dest='dmean', type=float, default=2.5, help='index of dispersion mean (default=2.5)')
parser.add_argument('--dispersion_sd', dest='dsd', type=float, default=0, help='index of dispersion sd (default=0)')
parser.add_argument('--seed', dest='seed', type=int, default=0, help='seed for random number generator')
args = parser.parse_args()

numpy.random.seed(args.seed)
epsilon = args.epsilon
def calculate_pl(dp_ref, dp_alt):
    """
    Calculate PL field in vcf file from allelic DP
    """
    dp = dp_ref + dp_alt
    if dp > 0:
        lik_homref = math.log(1 - epsilon, 10) * dp_ref + math.log(epsilon / 3, 10) * dp_alt
        lik_het = math.log((3 - 2 * epsilon) / 6.0, 10) * dp
        lik_homnonref = math.log(1 - epsilon, 10) * dp_alt + math.log(epsilon / 3, 10) * dp_ref
        lik = [lik_homref, lik_het, lik_homnonref]
        maxlik = max(lik)
        pl = [int(abs(10 * (l - maxlik))) for l in lik]
    else:
        pl = None
    return(pl)


rand = random.Random(0)
lamb = args.lamb
for line in sys.stdin:
    d = rand.normalvariate(mu=args.dmean, sigma=args.dsd)
    if d < 1:
        d = 1.001
    line = line.strip()
    if line[:2] == '##':
        print line
        continue
    elif line[:6] == '#CHROM':
        print '##analysis=simulate_dp.py --lambda %f --epsilon %f --dispersion_mean %f --dispersion_sd %f --seed %s' % (args.lamb, args.epsilon, args.dmean, args.dsd, str(args.seed))
    fields = line.split()
    genotypes = fields[9:]
    samples = []
    for gt in genotypes:
        # scipy nbinom takes (n, p) as arguments.
        # Convolution: sum_{i=1}{x}nbinom(n, p) = nbinom(xn, p).
        if gt.count('0')==0:
            dp_ref = 0
        else:
            dp_ref = nbinom.rvs(lamb * gt.count('0') / ((d - 1) * 2), 1 / d)
        if gt.count('1')==0:
            dp_alt = 0
        else:
            dp_alt = nbinom.rvs(lamb * gt.count('1') / ((d - 1) * 2), 1 / d)
        dp = dp_ref + dp_alt
        pl = calculate_pl(dp_ref, dp_alt)
        if pl:
            sample = "%s:%i,%i:%i:%i,%i,%i" % (gt, dp_ref, dp_alt, dp, pl[0], pl[1], pl[2])
            samples.append(sample)
        else:
            samples.append(gt)
    info = '%s;Dispersion=%s' % (fields[7], '{0:.3f}'.format(d))
    output = fields[:7] + [info] + ['GT:AD:DP:PL'] + samples
    print '\t'.join(output)
