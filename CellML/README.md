# AtrialModel_2025_Human/CellML
CellML and Python code for a novel muscle model based on non-diabetic and diabetic human data. This code has been tested on OpenCOR 0.8.2 and should run on at least versions from 0.7.1 onward. 

## Musgrave et al. (2025) Human atrial muscle model code

'musgrave_2025_diabetic.cellml' and 'musgrave_2025_nondiabetic.cellml' are the core model code that can be run directly to simulate isometric twitches. These each have a respective .sedml file which encodes the correct step size for the simulation as well as some standard plots.

## Jupyter Notebook

Opening an OpenCOR Jupyter kernel in this directory will allow the Human Cross-Bridge Modelling.ipynb notebook to be run (Instructions available in the notebook and [here](https://models.physiomeproject.org/e/afd). This notebook performs the basic simulations with the CellML models, but also allows implementation of work-loop simulations. These call on 3 additional modes of the baseline models, '..._isoton.cellml', '..._RS1.cellml' and '..._RS2.cellml', to implement the different types of length control required at different phases of the work-loop simulation

## Calcium input

All the CellML files use a measured calcium transient to drive the model. These transients are encoded in 'diabetic_calcium.cellml' and 'non-diabetic_calcium.cellml' and imported into each of the model files.

The raw data is included in csv form ('d_ca.csv' and 'nd_ca.csv'), as well as a simple Python script ('piecewise_Ca.py') that can be used to convert any array of time and calcium data into a piecewise description suitable for direct use in CellML code.