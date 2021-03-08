import numpy as np
import emcee
from multiprocessing import Pool
from NISSEN_4 import mw
import _pickle as cPickle


## Load the necessary data
#def load_data():
#    higha = np.genfromtxt('HIGH_no_binary.txt')
#    nhigh = len(higha[:, 0])

#    lowa = np.genfromtxt('LOW.txt')
#    nlow = len(lowa[:, 0])

#    ndata = nhigh + nlow
#    data = np.zeros((ndata, 6))
    
#    data[0:nhigh, 0] = higha[:, 10]
#    data[0:nhigh, 1] = higha[:, 12]
#    data[0:nhigh, 2] = higha[:, 15]
#    data[0:nhigh, 3] = higha[:, 16]
#   data[0:nhigh, 4] = higha[:, 8]/1000.
#    data[0:nhigh, 5] = higha[:, 9]/1000.
    
#    data[nhigh:ndata, 0] = lowa[:, 10]
#    data[nhigh:ndata, 1] = lowa[:, 12]
#    data[nhigh:ndata, 2] = lowa[:, 15]
#    data[nhigh:ndata, 3] = lowa[:, 16]
#    data[nhigh:ndata, 4] = lowa[:, 8]/1000.
#    data[nhigh:ndata, 5] = lowa[:, 9]/1000.
#    return data


# Logarithm of the prior
def lnprior(theta):
    at1, at2, atauin, aSDR = theta

    if (0.0 < at1 < 1.0 and 0.0 < at2 < 28.0 and 0.0 < atauin < 14.0 and 
         1.0 < aSDR < 50.):
        return 0.0

    return -np.inf


# Logarithm of the likelihood
def lnlikelihood(theta):
    at1, at2, atauin, aSDR = theta
    attspi, aafehv, aaspiv, aamfv,aAOFEv, aAMGFEv = mw(at1, at2, atauin, aSDR)

    lnlike = 0.
    for i in range(len(feh)):

        s = np.sqrt( ((feh[i] - aafehv[attspi>1.e-14])**2)/(efeh[i]**2) +
                          ((OFE[i] - aAOFEv[attspi>1.e-14])**2)/(eOFE[i]**2)+
                          ((MGFE[i] - aAMGFEv[attspi>1.e-14])**2)/(eMGFE[i]**2)+
                          ((age[i] - (13.7-attspi[attspi>1.e-14]))**2)/(eage[i]**2))
       
        ind = np.argmin(s)
        lnlike -= np.log((2. * np.pi)**(4/2) * efeh[i] * eOFE[i] *  eMGFE[i] *   eage[i])
        lnlike -= 0.5 * np.power((feh[i] - aafehv[ind]) / efeh[i], 2)
        lnlike -= 0.5 * np.power((OFE[i] - aAOFEv[ind]) / eOFE[i], 2)
        lnlike -= 0.5 * np.power((MGFE[i] - aAMGFEv[ind]) / eMGFE[i], 2)
        lnlike -= 0.5 * np.power((age[i] - (13.7 - attspi[ind])) / eage[i], 2)

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
ndim, nwalkers, nsteps = 4, 100, 1000



 


# Initial guess for the parameters
theta0 = (0.05, 10.0, 5.,  3.)


# Load the necessary data
#data = load_data()
global feh, efeh, OFE, eOFE, MGFE, eMGFE,age, eage

data = np.genfromtxt('NISSEN.dat')
feh, efeh = data[:, 3], data[:, 4]
OFE, eOFE = data[:, 5], data[:, 6]
MGFE, eMGFE = data[:, 7], data[:, 8]
age, eage = data[:, 1], data[:, 2]


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
