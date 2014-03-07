cdef extern from "math.h":
    double log (double x)
    double exp (double x)

from scipy.misc import comb
#from math import log
#from math import exp
from collections import namedtuple

# Defined as ('A' count, 'a', count, '-' count) (see GBStools notes).
GENO = ((2,0,0),
        (1,0,1),
        (1,1,0),
        (0,2,0),
        (0,1,1),
        (0,0,2))

# Make a hash of combinations to avoid redundant calculation.
choose = {}
for i in range(300):
    choose[i] = {}
    for j in range(i+1):
        choose[i][j] = comb(i, j, exact=1)

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

def trio_genotypes(genotypes, loglik=True):
    '''Calculate the probability of offspring genotypes in a trio.'''
    prob = {}
    for maternal_g in genotypes:
        for paternal_g in genotypes:
            maternal_gametes = []
            paternal_gametes = []
            gt = ParentalGT(maternal_g, paternal_g)
            prob[gt] = {}
            maternal_gametes = gametes(maternal_g)
            paternal_gametes = gametes(paternal_g)
            for m in maternal_gametes:
                for p in paternal_gametes:
                    offspring_g = tuple([a + b for a, b in zip(m, p)])
                    try:
                        prob[gt][offspring_g] += 0.25
                    except:
                        prob[gt][offspring_g] = 0.25
    if loglik:
        for gt in prob:
            for offspring_gt in prob[gt]:
                prob[gt][offspring_gt] = log(prob[gt][offspring_gt])
    return(prob)

# Get all possible parent/offspring combinations and their probabilities.
trio_gt = trio_genotypes(GENO)

def update(param, calls, disp):
    '''Update parameters by EM.'''
    n = len(calls)    # Number of samples.
    # Initialize counter variables for phi, delta, lambda, loglik.
    phi = [0, 0, 0]
    delta = 0
    lamb_numer = 0
    lamb_denom = 0
    loglik = 0
    exp_phi = {}
    exp_delta = {}
    for call in calls:
        sample_lik = {}
        sample_psi = {}
        exp_phi[call.sample] = [0, 0, 0]
        exp_delta[call.sample] = 0
        for z in (0, 1):
            for g in GENO:
                m = (2.0 - g[2])    # Apparent ploidy, given g.
                mu = param['lambda'] * call.NF * z * m / 2    # Expected coverage.
                psi = mu / (disp - 1)    # Size parameter for nbinom.
                D_lik = dreads(g, call.PL)    # Read data likelihood.
                d_lik = dnbinom(call.DP, mu, psi)    # Coverage likelihood.
                g_lik = dmultinom(g, param['phi'])    # Genotype likelihood.
                z_lik = dbernoulli(z, param['delta'])    # Digest fail likelihood.
                sample_lik[(z, g)] = D_lik + d_lik + g_lik + z_lik
                sample_psi[(z, g)] = psi
        # Normalize the likelihoods by likmax to avoid underflow.
        likmax = max(sample_lik.values())
        liksum = sum([exp(l - likmax) for l in sample_lik.values()])
        for z, g in sample_lik:
            psi = sample_psi[(z, g)]
            normlik = exp(sample_lik[(z, g)] - likmax)
            # Update the counter variables.
            phi[0] += g[0] * normlik / (2 * n * liksum)
            phi[1] += g[1] * normlik / (2 * n * liksum)
            phi[2] += g[2] * normlik / (2 * n * liksum)
            delta += (1 - z) * normlik  / (n * liksum)
            lamb_numer += (lambda_numer(g, z, call.DP, call.NF, param['lambda'], psi) *
                           normlik / liksum)
            lamb_denom += (lambda_denom(g, z, call.DP, call.NF, param['lambda'], psi) *
                           normlik / liksum)
            exp_phi[call.sample][0] += g[0] * normlik / liksum
            exp_phi[call.sample][1] += g[1] * normlik / liksum
            exp_phi[call.sample][2] += g[2] * normlik / liksum
            exp_delta[call.sample] += (1 - z) * normlik / liksum
        try:
            loglik += log(sum([exp(l) for l in sample_lik.values()]))
        except:
            loglik = None
    lamb = param['lambda'] - lamb_numer / lamb_denom
    if lamb < 0:
        lamb = min(1.0, param['lambda'] / 2)
    param_update = {'phi':phi, 
                    'delta':delta,
                    'lambda':lamb, 
                    'fail':param['fail'], 
                    'loglik':loglik,
                    'exp_phi':exp_phi,
                    'exp_delta':exp_delta}
    return(param_update)


def ped_update(param, calls, disp, parental_gt):
    '''Update parameters by EM for nuclear family.'''
    # Get hash of F1 GT prob, keyed by GT.
    offspring_gt = trio_gt[parental_gt]
    n = len(calls)
    delta = 0
    lamb_numer = 0
    lamb_denom = 0
    loglik = 0
    for call in calls:
        sample_lik = {}
        for z in (0, 1):
            if call.is_mother:
                genotypes = {parental_gt.mother:log(1)}
            elif call.is_father:
                genotypes = {parental_gt.father:log(1)}
            elif call.is_child:
                genotypes = offspring_gt
            else:
                genotypes = {}
            for g in genotypes:
                m = (2.0 - g[2])    # Apparent ploidy, given g.
                mu = param['lambda'] * call.NF * z * m / 2    # Expected coverage.
                psi = mu / (disp - 1)    # Size parameter for nbinom.
                D_lik = dreads(g, call.PL)    # Read data likelihood.
                d_lik = dnbinom(call.DP, mu, psi)    # Coverage likelihood.
                g_lik = genotypes[g]    # Genotype likelihood.
                z_lik = dbernoulli(z, param['delta'])    # Digest fail likelihood.
                sample_lik[(z, g, psi)] = D_lik + d_lik + g_lik + z_lik
        # Normalize the likelihoods by likmax to avoid underflow.
        likmax = max(sample_lik.values())
        if likmax == -float('Inf'):
            likmax = 0
            liksum = 1
        else:
            liksum = sum([exp(l - likmax) for l in sample_lik.values()])
        for z, g, psi in sample_lik:
            lik = exp(sample_lik[(z, g, psi)] - likmax)
            delta += (1 - z) * lik / (n * liksum)
            lamb_numer += (lambda_numer(g, z, call.DP, call.NF, param['lambda'], psi) *
                           lik / liksum)
            lamb_denom += (lambda_denom(g, z, call.DP, call.NF, param['lambda'], psi) *
                           lik / liksum)
        try:
            loglik += log(sum([exp(l) for l in sample_lik.values()]))
        except:
            loglik = -float('Inf')
    try:
        lamb = param['lambda'] - lamb_numer / lamb_denom
    except:
        lamb = param['lambda']
    param_update = {'delta':delta, 
                    'lambda':lamb, 
                    'fail':param['fail'], 
                    'loglik':loglik}
    return(param_update)     


cdef dreads(g, pl, loglik=True):
    # Function to calculate likelihood of observing reads given an apparent genotype.
    if g[2] == 0:
        g_apparent = g[1]
    elif g[2] == 1:
        g_apparent = 2 * g[1]
    else:
        g_apparent = None
    if pl:
        if g_apparent is not None:
            p = 10**(-float(pl[g_apparent]) / 10)
        else:
            p = 0.0
    else:
        p = 1.0
    if not loglik:
        return(p)
    else:
        if p > 0:
            return(log(p))
        else:
            return(-float('Inf'))

cdef dnbinom(x, mu, psi, loglik=True):
    '''Function to calculate negative binomial probability.'''
    psi = int(psi)    # Round to nearest integer to speed up comb().
    if psi <= 0:    # psi must be > 0.
        psi = 1
    if x:
        try:
            p = (choose[x + psi - 1][psi - 1] * 
                 (psi/(float(mu) + psi))**psi * 
                 (float(mu)/(float(mu) + psi))**x)
        except:
            p = (comb(x + psi - 1, psi - 1, exact=True) * 
                 (psi/(float(mu) + psi))**psi * 
                 (float(mu)/(float(mu) + psi))**x)
    else:
        p = (psi/(float(mu) + psi))**psi
    if not loglik:
        return(p)
    else:
        if p > 0:
            return(log(p))
        else:
            return(-float('Inf'))


cdef dmultinom(x, prob, loglik=True):
    '''Function to calculate multinomial probability.'''
    p = 2**(1 in x) * prob[0]**x[0] * prob[1]**x[1] * prob[2]**x[2]
    if not loglik:
        return(p)
    else:
        if p > 0:
            return(log(p))
        else:
            return(-float('Inf'))

cdef dbernoulli(x, prob, loglik=True):
    '''Function to calculate bernoulli probability.'''
    p = (1 - prob)**x * prob**(1 - x)
    if not loglik:
        return(p)
    else:
        if p > 0:
            return(log(p))
        else:
            return(-float('Inf'))

cdef lambda_numer(g, z, d, r, lamb, psi):
    '''Function to calculate numerator in EM lambda update (see GBStools notes).'''
    m = 2.0 - g[2]
    if m * z == 0:
        try:
            val = d / lamb
        except:
            val = 0
    else:
        try:
            val = d / lamb - (d + psi) / (lamb + 2 * psi / (r * z * m))
        except:
            val = 0
    return(val)


cdef lambda_denom(g, z, d, r, lamb, psi):
    '''Function to calculate denominator in EM lambda update (see GBStools notes).'''
    m = 2.0 - g[2]
    if m * z == 0:
        try:
            val = - d / lamb**2
        except:
            val = 0
    else:
        try:
            val = ((r * z * m / 2)**2 * (d + psi) / (lamb * r * z * m / 2 + psi)**2 -
                   d / lamb**2)
        except:
            val = 0
    return(val)
