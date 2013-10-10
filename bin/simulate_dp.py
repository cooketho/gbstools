import random
import sys
from scipy.stats import nbinom
import math
import numpy
from optparse import OptionParser

"""Use GT data in a GBS vcf file to simulate DP, PL, and AD data.
"""
USAGE = """
cat <myvcf> | simulate_geno.py -d <coverage> > <vcf with dp>
"""

parser = OptionParser(USAGE)
parser.add_option('-l', '--lambda', dest='lamb', type="int", help='mean depth of coverage')
(opt, args) = parser.parse_args()

numpy.random.seed(0)
epsilon = 0.001
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
d = 2.5
lamb = float(opt.lamb)
for line in sys.stdin:
    line = line.strip()
    if line[0] == '#':
        print line
        continue
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
    output = fields[:8] + ['GT:AD:DP:PL'] + samples
    print '\t'.join(output)
