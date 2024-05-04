# LS_WASP_IRT
This repository contains R codes for the paper titled “Optimizing Large-Scale Educational Assessment with a ‘Divide-and-Conquer’ Strategy: Fast and Efficient Distributed Bayesian Inference in IRT Models”.

The “2PL_model” folder and the “M2PL_model” folder contain the code for the location-scatter Wasserstein posterior (LS-WASP) algorithm and the PG-Gibbs (Polya-Gamma) algorithm in Section 3, respectively. Each file in the “2PL_model” folder include:
-The “Full_data_2PL(PG_Gibbs).R” file is the PG-Gibbs (Polya-Gamma) algorithm for the full data and the corresponding simulations for the 2PL model.
-The “subsets_2PL(LS_WASP).R” file is the LS-WASP algorithm for the subsets data and the corresponding simulations for the 2PL model. 

All codes are runnable and can be used to generate simulation results in subsection 5.1. The “M2PL_model” folder is structured similarly to the “2PL_model” folder. The code in the “M2PL_model” folder can be used to generate simulation results for subsection 5.2.
