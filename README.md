# Two-infall chemical evolution models
This repository can be used to fit the two-infall chemical evolution models to the observed stellar abundances and ages in a Bayesian framework based on Markov Chain Monte Carlo methods. For further details, see the following paper by Spitoni et al. (2020).
https://www.aanda.org/articles/aa/pdf/2020/03/aa37275-19.pdf

**Prerequisites**  
1. Usual python libraries  
1. *emcee*: The MCMC Hammer [http://dfm.io/emcee/current/]  
1. *corner* [https://corner.readthedocs.io/en/latest/]

**Troubleshooting**  
You may need to compile *MW_exp_4par.f*. Run the following in the command prompt:  
*f2py -c -m MW_exp_4par MW_exp_4par.f*

**Instructions**  
Check *MW_exp_4par.py* for the details of the fitting procedure. You may try running this script, however with default settings it may take upto a week on a modern computer with 8 threads. For quick tests, you may reduce the number of walkers and steps.
