# evolution-of-milky-way

**Prerequisites**  
1. Usual python libraries  
1. *emcee*: The MCMC Hammer [http://dfm.io/emcee/current/]  
1. *corner* [https://corner.readthedocs.io/en/latest/]

**Troubleshooting**  
You may need to compile *evolution.f*. Run the following in the command prompt:  
*f2py -c -m evolution evolution.f*

**Instructions**  
Look into the file *fitting.py* for the details of fitting.
