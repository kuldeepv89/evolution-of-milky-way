import numpy as np
import emcee
from multiprocessing import Pool
from MW_exp_3par_spitoni2019 import mw
import _pickle as cPickle


# Load the necessary data
def load_data():
    higha = np.genfromtxt('HIGH_no_binary.txt')
    nhigh = len(higha[:, 0])

    lowa = np.genfromtxt('LOW.txt')
    nlow = len(lowa[:, 0])

    ndata = nhigh + nlow
    data = np.zeros((ndata, 4))
    
    data[0:nhigh, 0] = higha[:, 10]
    data[0:nhigh, 1] = higha[:, 12]
    data[0:nhigh, 2] = higha[:, 15]
    data[0:nhigh, 3] = higha[:, 16]
#    data[0:nhigh, 4] = higha[:, 8]/1000.
#    data[0:nhigh, 5] = higha[:, 9]/1000.
    
    data[nhigh:ndata, 0] = lowa[:, 10]
    data[nhigh:ndata, 1] = lowa[:, 12]
    data[nhigh:ndata, 2] = lowa[:, 15]
    data[nhigh:ndata, 3] = lowa[:, 16]
#    data[nhigh:ndata, 4] = lowa[:, 8]/1000.
#    data[nhigh:ndata, 5] = lowa[:, 9]/1000.
    return data


# Logarithm of the prior
def lnprior(theta):
    at1, at2, atauin = theta

    if (0.0 < at1 < 7.0 and 0.0 < at2 < 14.0 and 0.0 < atauin < 14.0):
        return 0.0

    return -np.inf


# Logarithm of the likelihood
def lnlikelihood(theta):
    at1, at2, atauin = theta
    attspi, aafehv, aaspiv, aamfv = mw(at1, at2, atauin)

    lnlike = 0.
    for i in range(len(feh)):
        s = np.sqrt(np.power((feh[i] - aafehv[attspi>1.e-14]) / efeh[i], 2) + 
                    np.power((alh[i] - aaspiv[attspi>1.e-14]) / ealh[i], 2))
        ind = np.argmin(s)
        lnlike -= np.log((2. * np.pi) * efeh[i] * ealh[i])
        lnlike -= 0.5 * np.power((feh[i] - aafehv[ind]) / efeh[i], 2)
        lnlike -= 0.5 * np.power((alh[i] - aaspiv[ind]) / ealh[i], 2)
      

    return lnlike 


# Logarithm of the posterior probability
def lnprob(theta):
    lp = lnprior(theta)

    if not np.isfinite(lp):
        return -np.inf

    return lp + lnlikelihood(theta)



########### Start of the main program ################
# ndim     : number of free parameters
# nwalkers : number of walkers
# nsteps   : number of steps
ndim, nwalkers, nsteps = 3, 100, 1000


# Initial guess for the parameters
theta0 = (0.05,  5.,  3.)


# Load the necessary data
data = load_data()
global feh, efeh, alh, ealh
feh, efeh = data[:, 0], data[:, 1]
alh, ealh = data[:, 2], data[:, 3]
#age, eage = data[:, 4], data[:, 5]


# Parallelize the emcee run
with Pool() as pool:

    # Initialize the position of walkers
    pos = [theta0 + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]


    # Run MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
    sampler.run_mcmc(pos, nsteps, progress=False)


    # Print the diagnostics of the MCMC run
    # 0.2 < Mean acceptance fraction < 0.5
    # Mean autocorrelation time < nsteps
    print(('Mean acceptance fraction: %.3f') %(np.mean(sampler.acceptance_fraction)))
    #print(('Mean autocorrelation time: %.3f') %(np.mean(sampler.get_autocorr_time(tol=0))))


    # Write the chain
    with open('chain.pkl', 'wb') as fp:
        cPickle.dump(sampler.chain, fp)
