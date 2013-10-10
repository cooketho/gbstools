from optparse import OptionParser
import random
import numpy
from collections import namedtuple

"""Draw random genotypes for two parents, M and F, including
dropout alleles, and draw random genotypes for their offspring.
Output in VCF format.
"""

# Parse command line arguments from user.
USAGE = """
simulate_ped_vcf.py --samples <number of samples>
                    --sites <number of sites>
"""

parser = OptionParser(USAGE)
parser.add_option('-n', '--samples', dest='n', type='int', help='number of samples')
parser.add_option('-N', '--sites', dest='N', type='int', help='number of sites')
(opt, args) = parser.parse_args()

n = opt.n
N = opt.N
rand = random.Random(0)
numpy.random.seed(0)

# Defined as ('A' count, 'a', count, '-' count) (see GBStools notes).
GENOTYPES = ((2,0,0),
             (1,0,1),
             (1,1,0),
             (0,2,0),
             (0,1,1),
             (0,0,2))

GT = {(2,0,0):'0/0',
      (1,0,1):'0/.',
      (1,1,0):'0/1',
      (0,2,0):'1/1',
      (0,1,1):'1/.',
      (0,0,2):'./.'}

def gametes(g):
    '''Return a list of the possible gametes for genotype g.'''
    gametes = []
    for i in range(len(g)):
        gamete = [0] * len(g)
        gamete[i] = 1
        if g[i] == 1:
            gametes.append(gamete)
        elif g[i] == 2:
            gametes.append(gamete)
            gametes.append(gamete)
    return(gametes)

ParentalGT = namedtuple('ParentalGT', 'mother, father')

def trio_genotypes(g):
    '''Calculate the probability of offspring genotypes in an F1 cross.'''
    prob = {}
    for maternal_g in g:
        for paternal_g in g:
            maternal_gametes = []
            paternal_gametes = []
            gt = ParentalGT(maternal_g, paternal_g)
            prob[gt] = {}
            maternal_gametes = gametes(maternal_g)
            paternal_gametes = gametes(paternal_g)
            for m in maternal_gametes:
                for p in paternal_gametes:
                    f1_g = tuple([a + b for a, b in zip(m, p)])
                    try:
                        prob[gt][f1_g] += 0.25
                    except:
                        prob[gt][f1_g] = 0.25
    return(prob)

TRIOGENOTYPES = trio_genotypes(GENOTYPES)

header = "##fileformat=VCFv4.1\n"
header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
header += '\t'.join(['M', 'F'] + [str(i) for i in range(n)])
print header

sitenum = 0
snp = 0
for i in range(N):
    sitenum += 1
    parents = rand.choice(TRIOGENOTYPES.keys())
    if (0,0,2) not in parents:
        f1_gt = TRIOGENOTYPES[parents].keys()
        f1_prob = TRIOGENOTYPES[parents].values()
        f1_counts = numpy.random.multinomial(n, f1_prob)
        genotypes = []
        for gt, count in zip(f1_gt, f1_counts):
            genotypes += [gt] * count
        random.shuffle(genotypes)
        genotypes = [parents.mother, parents.father] + genotypes
        gt = [GT[i] for i in genotypes]
        dcount = sum([g[2] for g in genotypes])
        ac_sampled = sum([g[1] for g in genotypes])
        missing = gt.count('./.')
        if missing <= 0.2 * n:
            snp += 1
            info = "DCount=%i;ACgbs=%i;Missing=%i" % (dcount, ac_sampled, missing)
            fields = ['1', str(sitenum), str(snp)] + ['.'] * 4 + [info, 'GT'] + gt
            print '\t'.join(fields)
