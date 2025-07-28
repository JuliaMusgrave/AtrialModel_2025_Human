# AtrialModel_2025_Human
Simulation code for a novel cross-bridge and muscle model based on non-diabetic and diabetic human data, presented in [Musgrave et al., 2025](http://doi.org/10.1113/JP288463). Code has been provided in two languages. Note that a README is included within each folder to describe the contents in more detail.

## MATLAB

This is the original model code, containing both the cross-bridge version and the full muscle model version. This folder contains one script to run all of the figures in the publication and a second script to demonstrate use of the models more generally. mat files contain the different parameters needed to run non-diabetic and diabetic versions of the model.

## CellML

CellML code for each of the non-diabetic and diabetic muscle models is provided. This can be run directly to simulate isometric twitches. A Jupyter Notebook, to be run through OpenCOR, is also included to enable simulations of work-loops using several different modes for this model.

