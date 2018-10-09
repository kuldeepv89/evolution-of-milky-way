import numpy as np
import emcee
import corner
from evolution import mw
import _pickle as cPickle
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# Load the data
def load_data():

    higha = np.genfromtxt('HIGH.txt')
    nhigh = len(higha[:, 0])

    lowa = np.genfromtxt('LOW.txt')
    nlow = len(lowa[:, 0])

    ndata = nhigh + nlow
    data = np.zeros((ndata, 4))
    data[0:nhigh, 0] = higha[:, 10]
    data[0:nhigh, 1] = higha[:, 12]
    data[0:nhigh, 2] = higha[:, 15]
    data[0:nhigh, 3] = higha[:, 16]
    data[nhigh:ndata, 0] = lowa[:, 10]
    data[nhigh:ndata, 1] = lowa[:, 12]
    data[nhigh:ndata, 2] = lowa[:, 15]
    data[nhigh:ndata, 3] = lowa[:, 16]

    return data


# Logarithm of prior
def lnprior(theta):

    at1, at2, atauin = theta
    if 0.05 < at1 < 3.0 and 4.0 < at2 < 10.0 and 0.0 < atauin < 6.0:
        return 0.0

    return -np.inf


# Logarithm  of likelihood
def lnlikelihood(theta, feh, efeh, alh, ealh):

    at1, at2, atauin = theta
    attspi, aafehv, aaspiv, aamfv = mw(at1, at2, atauin)

    lnlike = 0.
    for i in range(len(feh)):
        s = np.sqrt((feh[i] - aafehv[attspi>1.e-14])**2 + (alh[i] - aaspiv[attspi>1.e-14])**2)
        ind = np.argmin(s)
        lnlike -= np.log(2. * np.pi * efeh[i] * ealh[i])
        lnlike -= 0.5 * np.power((feh[i] - aafehv[ind]) / efeh[i], 2)
        lnlike -= 0.5 * np.power((alh[i] - aaspiv[ind]) / ealh[i], 2)
        
    return lnlike 


# Logarithm of posterior probability
def lnprob(theta, feh, efeh, alh, ealh):

    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf

    return lp + lnlikelihood(theta, feh, efeh, alh, ealh)


# Sampling of posterior using emcee
def mcmc(ndim=3, nwalkers=100, nsteps=250, threads=6, theta0=(0.1, 8.0, 4.3)):

    data = load_data()
    feh, efeh, alh, ealh = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    theta = np.zeros(ndim)
    theta[0], theta[1], theta[2] = theta0[0], theta0[1], theta0[2]

    pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(feh, efeh, alh, ealh), threads=threads)

    for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
        print (("Iteration %d of %d completed...") %(i+1, nsteps))

    print(('Mean acceptance fraction: %.3f') %(np.mean(sampler.acceptance_fraction)))
    #print(('Autocorrelation time for at1, at2, atauin: %.3f %.3f %.3f') 
    #      %(sampler.acor[0], sampler.acor[1], sampler.acor[2]))

    with open('chain.pkl', 'wb') as fp:
        cPickle.dump(sampler.chain, fp)

    return


# Plot the time series using the emcee generated chain
def time_series(sampler_chain, fname='time_series.eps'):
    
    nwalkers, nsteps, ndim = sampler_chain.shape

    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Times New Roman')
    mpl.rc('text', usetex='false')
    mpl.rcParams.update({'font.size': 12})

    x = np.arange(nsteps)

    plt.figure(1)
    plt.subplot(3, 1, 1)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 0], '-', color='darkred')
    plt.ylabel(r'$\tau_{D1} \ ({\rm Gyr})$')
    plt.xlim(0, nsteps)
    plt.xticks([], [])
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))


    plt.subplot(3, 1, 2)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 1], '-', color='yellowgreen')
    plt.ylabel(r'$\tau_{D2} \ ({\rm Gyr})$')
    plt.xlim(0, nsteps)
    plt.xticks([], [])
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))

    plt.subplot(3, 1, 3)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 2], '-', color='darkblue')
    plt.xlabel('step number')
    plt.ylabel(r'$t_{\rm max} \ ({\rm Gyr})$')
    plt.xlim(0, nsteps)
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))

    plt.gcf().subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.10, hspace=0.05)
    plt.gcf().set_size_inches(6, 7)
    plt.savefig(fname)
    plt.close(1)

    return


# Make corner plot for fitting parameters using the emcee generated chain
def distribution(sampler_chain, nburnt=100, fname='corner.eps'):

    nwalkers, nsteps, ndim = sampler_chain.shape

    samples = sampler_chain[:, nburnt:, :].reshape((-1, ndim))

    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Times New Roman')
    mpl.rc('text', usetex='false')
    mpl.rcParams.update({'font.size': 12})

    plt.figure(1)
    fig = corner.corner(samples, color='darkblue', labels=[r'$\tau_{D1}$', r'$\tau_{D2}$', 
                        r'$t_{\rm max}$'], quantiles=[0.16, 0.5, 0.84], show_titles=True, 
                        title_kwargs={"fontsize": 12})
    fig.savefig(fname)
    plt.close(1)

    return


# Overplot evolution curves on [alpha/Fe] vs. [Fe/H] for all the walkers at final time step
def alphaFe_FeH(sampler_chain, fname='alphafe_feh.eps'):

    print ('Making [alpha/Fe] vs. [Fe/H] plot...')

    data = load_data()

    nwalkers, nsteps, ndim = sampler_chain.shape
    samples = sampler_chain[:, -1, :].reshape((-1, ndim))

    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Times New Roman')
    mpl.rc('text', usetex='false')
    mpl.rcParams.update({'font.size': 18})

    plt.figure(1)
    plt.errorbar(data[:, 0], data[:, 2], xerr=data[:, 1], yerr=data[:, 3], fmt='none', 
                 mfc='darkblue', mec='darkblue', capsize=0)
    plt.plot(data[:, 0], data[:, 2], 'o', ms=4, mfc='darkblue', mec='darkblue')
    for i in range(nwalkers):
        at1, at2, atauin = samples[i, 0], samples[i, 1], samples[i, 2]
        print (('Calculating track %d of %d...') %(i+1, nwalkers))
        attspi, aafehv, aaspiv, aamfv = mw(at1, at2, atauin)
        plt.plot(aafehv[attspi>1.e-14], aaspiv[attspi>1.e-14], '-', color='mistyrose')
    print ('Done...')        
    plt.xlabel(r'$[{\rm Fe}/{\rm H}] \ ({\rm dex})$')
    plt.ylabel(r'$[\alpha/{\rm Fe}] \ ({\rm dex})$')
    plt.xlim(-1.0, 0.5)
    plt.ylim(-0.1, 0.3)
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(7))

    plt.gcf().set_size_inches(8, 6)
    plt.savefig(fname)    
    plt.close(1)

    return




########### Start of main program ################
ndim, nwalkers, nsteps = 3, 10, 10


print (datetime.now())


########### Run the MCMC chain ###############
_ = mcmc(ndim=ndim, nwalkers=nwalkers, nsteps=nsteps, threads=6, theta0=(0.3, 8.0, 4.3))
print (datetime.now())


########### Load the stored chain #################
with open('chain.pkl', 'rb') as fp:
    sampler_chain = cPickle.load(fp) 


########### Make the time series plot ##################
_ = time_series(sampler_chain, fname='time_series.eps')


########### Make the corner plot #################
_ = distribution(sampler_chain, nburnt=0, fname='corner.eps')


########### Make [Fe/H] vs. [alpha/H] plot ################
_ = alphaFe_FeH(sampler_chain, fname='alphafe_feh.eps')


print (datetime.now())
