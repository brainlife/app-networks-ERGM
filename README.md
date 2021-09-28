# app-networks-ERGM
Exponential Random Graph Models (ERGMs) for brain networks estimated via spectral coherence between high-density electroencephalographic (EEG) signals

 The following are the inputs outputs :
IN - data/brainM.txt : a (NxN) adiacency matrix, representing measured correlations.
IN - data/attrs.txt [optional] : a (Nxm) matrix, each column represents an attribute for the N nodes, m is the total number of attributes.
IN - data/ecov_1.txt [optional] : a (NxN) matrix representing a edge covariate for the N^2 possible edges.

OUT - output/mcmc-diagnostic.pdf : a pdf with some plots
OUT - output/gof.pdf : a pdf with some plots
OUT - output/log_computation.txt : log of the computation sunk from the terminal
OUT - output/estimation.txt : print of the estimated parameters sunk from the terminal
