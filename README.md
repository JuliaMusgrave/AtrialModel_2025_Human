# AtrialModel_2025_Human
MATLAB code for a novel cross-bridge and muscle model based on non-diabetic and diabetic human data. This code has been tested on MATLAB R2024a and should run on at least versions from 2021 onward.

## Musgrave et al. (2025) Human atrial muscle model code

Core model code is 'Mmodel_2025_Human.m', set up to be run with a MATLAB ode solver.

There is also a maximally-activated cross-bridge model 'XBmodel_2025_Human' set up similarly. 

## Scripts
Example script running the model (either non-diabetic or diabetic version) for isometric and work-loop simulations: 'run_human_Ca_transient.m'. This script contains the Ca50 parameters required for the non-diabetic and diabetic Ca transients provided. Also includes very basic use of the cross-bridge only model.

In addition, a script to reproduce the figures in the paper is provided, 'paper_figure_plotting.m'
This script shows the model fits to the average data ('average_human_fitting_data.mat', also available on [figshare](https://figshare.com/articles/dataset/Human_atrial_cross-bridge_data_diabetic_and_non-diabetic_/28180337?file=51622616)) and then performs a range of isometric and work-loop simulations for both diabetic and non-diabetic models. Also includes a sensitivity analysis. 

## Parameter data files
Several .mat files containing specific model parameters are provided and called in these scripts to run the model:

- ND_xb_fit/D_xb_fit - provide the cross-bridge model parameters based on non-diabetic and diabetic data. 'x_i' is the variable that takes into account the effect of the complete muscle model with Ca, while 'x_p' is based on the fully-activated cross-bridge model.

- ND_pass_fit/D_pass_fit - provide the passive model parameters based on non-diabetic and diabetic data

- thin_fil_ps - provides the thin filament/Ca regulation parameters, based on those used in the [Land 2017](https://pubmed.ncbi.nlm.nih.gov/28392437/) model. NOTE that the first parameter Ca50 is changed to account for differences in the Ca transient used to run the model

- Ca_transients - provides a description of the non-diabetic and diabetic Ca transients that can be used to run the model. These are provided in an array of points obtained at 1 kHz. Based on measurements reported in [Jones et al. 2023](https://physoc.onlinelibrary.wiley.com/doi/10.14814/phy2.15599)

## Other included functions
'XBmodel_2024_linear_perms.m' from [2024 rat model](https://github.com/JuliaMusgrave/XBModel_2024_Rat) is also included to run the linearised version of the XB model for fit visualisation. 

'full_sensitivity_for_paper.m' and 'WL_function.m' are helper functions to run sensitivity analysis and work-loop simulations, respectively. 'twitch_analysis.m' finds the amplitude and duration of a isometric twitch for analysis.

'get_sim_inputs.m' is a helper function for running the model and 'SSsim_par' is helper function for running the model until the force reaches steady-state.

