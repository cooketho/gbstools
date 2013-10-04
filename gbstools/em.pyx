cdef extern from "math.h":
    double log (double x)
    double exp (double x)

from scipy.misc import comb

# Make a hash of combinations to avoid redundant calculation.
choose = {}
for i in range(300):
    choose[i] = {}
    for j in range(i+1):
        choose[i][j] = comb(i, j, exact=1)

# Defined as ('A' count, 'a', count, '-' count) (see GBStools notes).
GENOTYPES1 = ((2,0,0),
              (1,0,1),
              (1,1,0),
              (0,2,0),
              (0,1,1),
              (0,0,2))
# For running EM without PL data.
GENOTYPES2 = ((0,0,0),
              (0,0,1),
              (0,0,2))

def update(param, calls, disp):
    '''Update parameters by EM.'''
    dp = [call.DP for call in calls]
    r = [call.NF for call in calls]
    # If there is no PL data, set pl to None.
    if [call.PL for call in calls if call.PL]:
        pl = [call.PL for call in calls]
    else:
        pl = None
    n = len(dp)    # Number of samples.
    if pl:
        genotypes = GENOTYPES1
    else:
        genotypes = GENOTYPES2
    # Initialize counter variables for phi, delta, lambda, loglik.
    phi = [0, 0, 0]
    delta = 0
    lamb_numer = 0
    lamb_denom = 0
    loglik = 0
    for i in range(n):
        sample_lik = {}
        for z in (0, 1):
            for g in genotypes:
                m = (2.0 - g[2])    # Apparent ploidy, given g.
                mu = param['lambda'] * r[i] * z * m / 2    # Expected coverage.
                psi = mu / (disp - 1)    # Size parameter for nbinom.
                if pl:
                    D_lik = dreads(g, pl[i])    # Read data likelihood.
                else:
                    D_lik = 0
                d_lik = dnbinom(dp[i], mu, psi)    # Coverage likelihood.
                g_lik = dmultinom(g, param['phi'])    # Genotype likelihood.
                z_lik = dbernoulli(z, param['delta'])    # Digest fail likelihood.
                sample_lik[(z, g, psi)] = D_lik + d_lik + g_lik + z_lik
        # Normalize the likelihoods by likmax to avoid underflow.
        likmax = max(sample_lik.values())
        liksum = sum([exp(l - likmax) for l in sample_lik.values()])
        for z, g, psi in sample_lik:
            lik = exp(sample_lik[(z, g, psi)] - likmax)
            # Update the counter variables.
            phi[0] += g[0] * lik / (2 * n * liksum)
            phi[1] += g[1] * lik / (2 * n * liksum)
            phi[2] += g[2] * lik / (2 * n * liksum)
            delta += (1 - z) * lik / (n * liksum)
            lamb_numer += (lambda_numer(g, z, dp[i], r[i], param['lambda'], psi) *
                           lik / liksum)
            lamb_denom += (lambda_denom(g, z, dp[i], r[i], param['lambda'], psi) *
                           lik / liksum)
        try:
            loglik += log(sum([exp(l) for l in sample_lik.values()]))
        except:
            loglik = None
    lamb = param['lambda'] - lamb_numer / lamb_denom
    param_update = {'phi':phi, 
                    'delta':delta, 
                    'lambda':lamb, 
                    'converged':param['converged'], 
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
    # Function to calculate negative binomial probability.
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
    # Function to calculate multinomial probability.
    p = 2**(1 in x) * prob[0]**x[0] * prob[1]**x[1] * prob[2]**x[2]
    if not loglik:
        return(p)
    else:
        if p > 0:
            return(log(p))
        else:
            return(-float('Inf'))

cdef dbernoulli(x, prob, loglik=True):
    # Function to calculate bernoulli probability.
    p = (1 - prob)**x * prob**(1 - x)
    if not loglik:
        return(p)
    else:
        if p > 0:
            return(log(p))
        else:
            return(-float('Inf'))

cdef lambda_numer(g, z, d, r, lamb, psi):
    # Function to calculate numerator in EM lambda update (see GBStools notes).
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
    # Function to calculate denominator in EM lambda update (see GBStools notes).
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
