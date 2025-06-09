# Analysis-and-modeling-of-voltage-gated-sodium-channels-clusters
The code used in the manuscript: Tarasov et al. Clustering Dynamically Modulate the Biophysics of Voltage-Gated Sodium Channels: How Nanoscale Phenomena Determine Health and Disease. bioRxiv 2025.05.31.657169; doi: https://doi.org/10.1101/2025.05.31.657169

**Cell-attached recordings analysis**: contains functions to (1) load patch-clamp recording from .abf file, (2) adjust baseline of current sweeps, (3) measure single-channel amplitudes by fitting gaussian mixture distributions to all current amplitudes, (4) subtracting capacitive transient currents, (5) de-noisation (idealization) of current sweeps using the Viterbi algorithm, (6) plotting current sweeps, (7) ensemble averaging idealized current sweeps, (8) calculation of mean late sodium current by normalization to peak sodium currents in idealized current sweeps, (9) fitting ensemble-average peak sodium current decay to the exponential functions to estimate the time constant of peak current decay.

**Quantification of functional interaction between channels**: contains functions to (1) fitting mean-variance relations in numbers of open channels to the binomial mean-variance relations, (2) calculation of Kullback–Leibler divergence and (3) variance-to-mean ratio (coefficient of dispersion) in numbers of open channels based on idealized current sweeps.

**Bayesian inference of the hidden Markov model parameters for currents idealization**: contains code for (1) loading, visualization and selection of single channel current sweeps used for subsequent inference, (2) Bayesian model for inference of the hidden Markov model parameters, (4) sampling the hidden Markov model parameters using tensorflow-probability’s NUTS sampler. 

**Single NaV1.5 channel Bayesian Moreno model optimization**: contains code for (1) Bayesian inference Moreno model parameters and (2) deterministic and stochastic simulations of Markov models.

**Composite models of NaV1.5 channels**: contains code for (1) building of composite models of non-interacting and interacting WT- and ΔKPQ-NaV1.5 channels in the absence and presence of lidocaine and (2) simulation of ion channel activity with these models.

**Random simulation of channel interaction in the membrane**: contains code for simulation of channel localizations in the membrane and prediction of integral sodium currents based on numbers of interacting and noninteracting channels in pairs given deterministic simulations of composite models of WT-NaV1.5.

**Stochastic simulation of Naundorf model of 2 channels**: contains code for stochastic simulation of activity of 2 interacting or non-interacting channels.

**Estimate of average single channel open probability in NaV1.5 clusters based on minflux and patch-clamp**: contains code to estimate differences in average single channel peak open probabilities in control and reduced surface expression conditions based on minflux measurements of cluster densities and ensemble average peak sodium currents in multi-channel cell-attached recordings.
