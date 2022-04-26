# bHUB
## Bayesian inference of hub nodes across multiple networks

Author: Junghi Kim

Contact: junghikim0@gmail.com



The given Matlab files for Bayesian inference of hub nodes for multiple networks (bHUB) are associated with the following publication:

Kim J, Do KA, Ha MJ and Peterson CB. (2019). Bayesian inference of hub nodes for multiple networks. *Biometrics*. 75(1): 172-182.


These scripts are related to the Matlab code for Stochastic Search Structure Learning (SSSL) provided by Hao Wang at https://msu.edu/~haowang/RESEARCH/Wang2013WP.pdf


Associated publication:

Wang H, Scaling It Up: Stochastic Search Structure Learning in Graphical Models. *Bayesian Analysis*. 10 (2015): 351-377. 


Please cite both publications if you use this code. Thanks!


## OVERVIEW OF FILES 

Start with the demo file run.m, which conducts the simulation experiments for an illustrative example in Kim et al. (2019)

## run.m
Basic example of running MCMC sampler and generating results summaries on a simple setting with
3 groups shown in an illustrative example


## hubGGMMCMC.m
Code for running MCMC sampler. Hyperparameter settings (which can be modified if needed)
are defined within this file


## hubGGM_SSVS.m
Helper function for updating precision matrix and graph using a modified version of SSSL


## calc_mrf_C.m
Helper function for calculating normalizing constant

Some additional helper functions are also provided.



