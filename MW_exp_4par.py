import numpy as np
import emcee
from multiprocessing import Pool
from MW_exp_4par import mw
import _pickle as cPickle
import corner
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# Load the necessary data
def load_data():
    higha = np.genfromtxt('HIGH_no_binary.txt')
    nhigh = len(higha[:, 0])

    lowa = np.genfromtxt('LOW.txt')
    lowa = lowa[lowa[:, 8]<=10000., :]
    nlow = len(lowa[:, 0])

    ndata = nhigh + nlow
    data = np.zeros((ndata, 6))
    
    data[0:nhigh, 0] = higha[:, 10]
    data[0:nhigh, 1] = higha[:, 12]
    data[0:nhigh, 2] = higha[:, 15]
    data[0:nhigh, 3] = higha[:, 16]
    data[0:nhigh, 4] = higha[:, 8]/1000.
    data[0:nhigh, 5] = higha[:, 9]/1000.
    
    data[nhigh:ndata, 0] = lowa[:, 10]
    data[nhigh:ndata, 1] = lowa[:, 12]
    data[nhigh:ndata, 2] = lowa[:, 15]
    data[nhigh:ndata, 3] = lowa[:, 16]
    data[nhigh:ndata, 4] = lowa[:, 8]/1000.
    data[nhigh:ndata, 5] = lowa[:, 9]/1000.
    return data


# Logarithm of the prior
def lnprior(theta):
    at1, at2, atauin, aSDR = theta

    if (0.0 < at1 < 7.0 and 0.0 < at2 < 28.0 and 0.0 < atauin < 14.0 and 
         1.0 < aSDR < 50.):
        return 0.0

    return -np.inf


# Logarithm of the likelihood
def lnlikelihood(theta):
    at1, at2, atauin, aSDR = theta
    attspi, aafehv, aaspiv, aamfv = mw(at1, at2, atauin, aSDR)

    lnlike = 0.
    for i in range(len(feh)):
        s = np.sqrt(np.power((feh[i] - aafehv[attspi>1.e-14]) / efeh[i], 2) +
                    np.power((alh[i] - aaspiv[attspi>1.e-14]) / ealh[i], 2) +
                    np.power((age[i] - (13.7 - attspi[attspi>1.e-14])) / eage[i], 2))
        ind = np.argmin(s)
        lnlike -= np.log((2. * np.pi)**(3/2) * efeh[i] * ealh[i] * eage[i])
        lnlike -= 0.5 * np.power((feh[i] - aafehv[ind]) / efeh[i], 2)
        lnlike -= 0.5 * np.power((alh[i] - aaspiv[ind]) / ealh[i], 2)
        lnlike -= 0.5 * np.power((age[i] - (13.7 - attspi[ind])) / eage[i], 2)

    return lnlike 


# Logarithm of the posterior probability
def lnprob(theta):
    lp = lnprior(theta)

    if not np.isfinite(lp):
        return -np.inf

    return lp + lnlikelihood(theta)


# Plot the time series using the emcee generated chain
def time_series(sampler_chain, fname='./time_series.eps'):
    
    nwalkers, nsteps, ndim = sampler_chain.shape

    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Times New Roman')
    mpl.rc('text', usetex='false')
    mpl.rcParams.update({'font.size': 12})

    x = np.arange(nsteps)

    plt.figure(1)
    plt.subplot(4, 1, 1)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 0], '-', color='darkred')
    plt.ylabel(r'$\tau_{D1} \ ({\rm Gyr})$')
    plt.xlim(0, nsteps)
    plt.xticks([], [])
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))


    plt.subplot(4, 1, 2)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 1], '-', color='darkgreen')
    plt.ylabel(r'$\tau_{D2} \ ({\rm Gyr})$')
    plt.xlim(0, nsteps)
    plt.xticks([], [])
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))

    plt.subplot(4, 1, 3)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 2], '-', color='darkblue')
    plt.xlabel('step number')
    plt.ylabel(r'$t_{\rm max} \ ({\rm Gyr})$')
    plt.xlim(0, nsteps)
    plt.xticks([], [])
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))

    plt.subplot(4, 1, 4)
    for i in range(nwalkers):
        plt.plot(x, sampler_chain[i, :, 3], '-', color='lime')
    plt.xlabel('step number')
    plt.ylabel(r'$\sigma_D/\sigma_T$')
    plt.xlim(0, nsteps)
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(5))


    plt.gcf().subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.10, hspace=0.05)
    plt.gcf().set_size_inches(6, 7)
    
    plt.savefig(fname)
    plt.close(1)

    return


# Make corner plot for fitting parameters using the emcee generated chain
def distribution(sampler_chain, nburn=100, fname='./corner.eps'):

    nwalkers, nsteps, ndim = sampler_chain.shape
    samples = sampler_chain[:, nburn:, :].reshape((-1, ndim))

    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Times New Roman')
    mpl.rc('text', usetex='false')
    mpl.rcParams.update({'font.size': 12})

    plt.figure(1)
    fig = corner.corner(samples, color='darkblue', labels=[r'$\tau_{D1}$', 
                        r'$\tau_{D2}$', r'$t_{\rm max}$', r'$\sigma_D/\sigma_T$'], 
                        quantiles=[0.16, 0.5, 0.84], show_titles=True, 
                        title_kwargs={"fontsize": 12})
    fig.savefig(fname)
    plt.close(fig)

    return


# Abundance and age plots 
def abundance_and_age_plots(sampler_chain, fname=['./alphafe_feh.eps', 
                            './alphafe_age.eps', './MH_age.eps']):

    data = load_data()

    nwalkers, nsteps, ndim = sampler_chain.shape
    samples = sampler_chain[:, -1, :].reshape((-1, ndim))

    attspi = np.zeros((3000, nwalkers))
    aafehv = np.zeros((3000, nwalkers))
    aaspiv = np.zeros((3000, nwalkers))
    aamfv = np.zeros((3000, nwalkers))
    for i in range(nwalkers):
        print (('Calculating track %d of %d...') %(i+1, nwalkers))
        at1, at2, atauin, aSDR  = (samples[i, 0], samples[i, 1], samples[i, 2], 
                                   samples[i, 3])
        attspi[:, i], aafehv[:, i], aaspiv[:, i], aamfv[:, i] = mw(at1, at2, atauin, aSDR)
    print ('Done...')        


    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Times New Roman')
    mpl.rc('text', usetex='false')
    mpl.rcParams.update({'font.size': 18})


    # Alpha vs. [Fe/H]
    plt.figure(1)
    plt.errorbar(data[:, 0], data[:, 2], xerr=data[:, 1], yerr=data[:, 3], fmt='none', 
                 mfc='darkblue', mec='darkblue', capsize=0)
    plt.plot(data[:, 0], data[:, 2], 'o', ms=4, mfc='darkblue', mec='darkblue')
    for i in range(nwalkers):
        plt.plot(aafehv[attspi[:, i]>1.e-14, i], aaspiv[attspi[:, i]>1.e-14, i], '-', 
                 color='mistyrose')
    plt.axhline(y=0., c='k', ls=':')
    plt.axvline(x=0., c='k', ls=':')
    plt.xlabel(r'$[{\rm Fe}/{\rm H}] \ ({\rm dex})$')
    plt.ylabel(r'$[\alpha/{\rm Fe}] \ ({\rm dex})$')
    plt.xlim(-1.0, 0.5)
    plt.ylim(-0.1, 0.3)
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(7))

    plt.gcf().set_size_inches(8, 6)
    plt.savefig(fname[0])    
    plt.close(1)


    # Alpha vs. age
    plt.figure(2)
    plt.errorbar(data[:, 4], data[:, 2], xerr=data[:, 5], yerr=data[:, 3], fmt='none', 
                 mfc='darkblue', mec='darkblue', capsize=0)
    plt.plot(data[:, 4], data[:, 2], 'o', ms=4, mfc='darkblue', mec='darkblue')
    for i in range(nwalkers):
        plt.plot(13.7 - attspi[attspi[:, i]>1.e-14, i], aaspiv[attspi[:, i]>1.e-14, i], 
                 '-', color='mistyrose')
    plt.xlabel(r'Age [Gyr]')
    plt.ylabel(r'$[\alpha/{\rm Fe}]$')
    plt.xlim(0, 15)
    plt.ylim(-0.05, 0.25)
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(7))

    plt.gcf().set_size_inches(8, 6)
    plt.savefig(fname[1])    
    plt.close(2)


    # [M/H] vs. age
    plt.figure(3)
    plt.errorbar(data[:, 4], data[:, 0], xerr=data[:, 5], yerr=data[:, 1], fmt='none', 
                 mfc='darkblue', mec='darkblue', capsize=0)
    plt.plot(data[:, 4], data[:, 0], 'o', ms=4, mfc='darkblue', mec='darkblue')
    for i in range(nwalkers):
        x = 13.7 - attspi[attspi[:, i]>1.e-14, i]
        y = (aafehv[attspi[:, i]>1.e-14, i] + 
             np.log10(0.638*10.**(aaspiv[attspi[:, i]>1.e-14, i]) + 0.362))
        plt.plot(x, y, '-', color='mistyrose')
    plt.xlabel(r'Age [Gyr]')
    plt.ylabel(r'$[{\rm M}/{\rm H}]$')
    plt.xlim(0, 15)
    plt.ylim(-1.5, 1)
    plt.gca().xaxis.set_major_locator(MaxNLocator(7))
    plt.gca().yaxis.set_major_locator(MaxNLocator(7))

    plt.gcf().set_size_inches(8, 6)
    plt.savefig(fname[2])    
    plt.close(3)

    return



########### Start of the main program ################
# ndim     : number of free parameters
# nwalkers : number of walkers
# nsteps   : number of steps
# nburn    : steps to burn
ndim, nwalkers, nsteps, nburn = 4, 100, 500, 400


# Initial guess for the parameters
theta0 = (0.05, 10.0, 5.,  3.)


# Load the necessary data
data = load_data()
global feh, efeh, alh, ealh, age, eage
feh, efeh = data[:, 0], data[:, 1]
alh, ealh = data[:, 2], data[:, 3]
age, eage = data[:, 4], data[:, 5]


# Parallelize the emcee run
with Pool() as pool:

    # Initialize the position of walkers
    pos = [theta0 + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]


    # Run MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
    sampler.run_mcmc(pos, nsteps, progress=True)


    # Print the diagnostics of the MCMC run
    # 0.2 < Mean acceptance fraction < 0.5
    # Mean autocorrelation time < nsteps
    print(('Mean acceptance fraction: %.3f') %(np.mean(sampler.acceptance_fraction)))
    #print(('Mean autocorrelation time: %.3f') %(np.mean(sampler.get_autocorr_time(tol=0))))


    # Write the chain
    with open('chain.pkl', 'wb') as fp:
        cPickle.dump(sampler.chain, fp)


# Load the stored chain
with open('chain.pkl', 'rb') as fp:
    sampler_chain = cPickle.load(fp)


# Make few useful plots
_ = time_series(sampler_chain)
_ = distribution(sampler_chain, nburn=nburn)
_ = abundance_and_age_plots(sampler_chain)
