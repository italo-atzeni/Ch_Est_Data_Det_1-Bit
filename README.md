# Channel Estimation and Data Detection Analysis of Massive MIMO with 1-Bit ADCs

This repository contains the MATLAB code to generate numerical results for the paper:

I. Atzeni and A. Tölli, "Channel Estimation and Data Detection Analysis of Massive MIMO with 1-Bit ADCs," *IEEE Trans. Wireless Commun. (to appear)*, 2021.

**Abstract of the paper.** We present an analytical framework for the channel estimation and the data detection in massive multiple-input multiple-output uplink systems with 1-bit analog-to-digital converters (ADCs) and i.i.d. Rayleigh fading. First, we provide closed-form expressions of the mean squared error (MSE) of the channel estimation considering the state-of-the-art linear minimum MSE estimator and the class of scaled least-squares estimators. For the data detection, we provide closed-form expressions of the expected value and the variance of the estimated symbols when maximum ratio combining is adopted, which can be exploited to efficiently implement minimum distance detection and, potentially, to design the set of transmit symbols. Our analytical findings explicitly depend on key system parameters such as the signal-to-noise ratio (SNR), the number of user equipments, and the pilot length, thus enabling a precise characterization of the performance of the channel estimation and the data detection with 1-bit ADCs. The proposed analysis highlights a fundamental SNR trade-off, according to which operating at the right noise level significantly enhances the system performance.

**Link to the paper on arXiv:** [https://arxiv.org/pdf/2102.10172.pdf](https://arxiv.org/pdf/2102.10172.pdf)

**Link to the paper on IEEEXplore:** TBD

# License and referencing

This code is licensed under the GPLv2 license. If you use this code in any way for research that results in publications, please cite the above article.

# Software and harware requirements

The code has been tested in MATLAB 2020b.

The numerical results based on the analytical functions can be generated smoothly (i.e., within seconds or minutes) on a MacBook Pro (processor: 2.4 GHz Dual-Core Intel Core i7, memory: 8 GB 1867 MHz LPDDR3). For the numerical results based on Monte Carlo simulations, we recommend using "parfor" on a computing server with at least 10 cores, especially for the file SER_VS_rho.m

# How to generate the figures

The main scripts to generate the figures are listed below. Each of these scripts saves a binary file in the "files_mat" folder. Note that the "files_mat" folder already contains all the binary files that are necessary for the plots.

All the other scripts are used to compute specific functions and variables defined throughout the paper.

**Figure 2**
- MSE_VS_rho.mat: Plot the MSE of the channel estimation against the transmit SNR.
- MSE_VS_rho_p1.mat: Plot the MSE of the channel estimation against the transmit SNR with K=1 and p=1_tau.

**Figure 3**
- MSE_VS_tau.mat: Plot the MSE of the channel estimation against the pilot length.
- MSE_VS_tau_opt_p1.mat: Plot the optimal MSE of the channel estimation and the optimal transmit SNR against the pilot length with K=1 and p=1_tau.

**Figure 4**
- dec_scatter.m: Plot the estimated symbols with the MRC receiver and the corresponding expected values.

**Figure 5**
- var_VS_M.m: Plot the normalized variance of the estimated symbols against the number of BS antennas.
- var_VS_M_p1.m: Plot the normalized variance of the estimated symbols against the number of BS antennas with p=1_tau.

**Figure 6**
- var_VS_rho.m: Plot the normalized variance of the estimated symbols against the transmit SNR.
- var_VS_rho_p1.m: Plot the normalized variance of the estimated symbols against the transmit SNR with p=1_tau.

**Figure 7**
- var_VS_tau.m: Plot the normalized variance of the estimated symbols against the pilot length.
- var_VS_tau_p2.m: Plot the normalized variance of the estimated symbols against the pilot length with p=1_tau.

**Figure 8**
- SER_VS_rho.m: Plot the SER against the transmit SNR.
