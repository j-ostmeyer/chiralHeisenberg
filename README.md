# Chiral Heisenberg Model in 2+1 dimensions with Domain Wall Fermions
## Simon Hands and Johann Ostmeyer

This repository contains the program required to reproduce the results presented in *"Lattice Field Theory Analysis of the Chiral Heisenberg Model"*, [PRB 114, 125102 (2026)](https://journals.aps.org/prb/abstract/10.1103/pndw-743y), [arXiv:2603.29461 [hep-lat]](https://arxiv.org/abs/2603.29461).

It contains a program-come-library for computations in lattice field theory with domain
wall fermions in 2+1D. Developed by Simon Hands.
Parallelised and features added by Johann Ostmeyer.

The data (raw as well as analysed) required to reproduce the results of the paper is published under open access [10.17638/datacat.liverpool.ac.uk/3076](https://doi.org/10.17638/datacat.liverpool.ac.uk/3076).

## Order Parameter Data
The raw data produced by the simulation program is found in the files
'condensates.csv'.  The order parameter |Phi| is defined in Eq. (13) of the
accompanying paper.  The data is arranged in rows for each ensemble generated.
The columns are labelled as follows:

L:  the linear extent of the spacetime volume L^3\
L_s:  the domain wall separation\
beta:  the inverse coupling parameter

no. meas: the number of measurements 

avg: mean order parameter |Phi|\
error:   standard error in |Phi|

tau_int:  integrated autocorrelation time\
acc:   RHMC acceptance rate

<phi^2>:   mean square order parameter\
d_<phi^2>:   standard error in <phi^2>

<(phi^2)^2>:  fourth moment of order parameter\
d_<(phi^2)^2>: standard error in <(phi^2)^2)>

bind. cum.:   Binder cumulant\
d_bind. cum: standard error in Binder cumulant

## Fermion Propagator Data


The correlators $S_\gamma$ and $S_m$ defined in Eq. (19) of the paper are written as complex numbers into the following output streams:
- $(S_\gamma)_{11}$ :	*fort.5011*
- $(S_\gamma)_{12}$ :	*fort.5012*
- $(S_m)_{11}$ :	*fort.5001*
- $(S_m)_{12}$ :	*fort.5002*

The streams are collected in directories named according to the ensemble parameters, e.g. `L_12_Ls_24_b_0.475`. All data is packed into the file 'fermion_data.tar'

Results from fits using the Truncated Hankel Correlator (THC) method are found in the files `thc.csv` and `thc_coeffs.csv`.
The columns are labelled as follows:

beta:  the inverse coupling parameter\
L,Lt:  the linear extent of the space-time volume $L^3=L^2 L_t$\
L_s:  the domain wall separation

channel: fit ansatz, only 6 (THC) relevant for this paper

eff. mass: ground state energy shown in Fig. 8 of the paper (`thc.csv`) or ground state contribution to the slope shown in Fig. 9 of the paper (`thc_coeffs.csv`)\
stat. err: corresponding statistical error
