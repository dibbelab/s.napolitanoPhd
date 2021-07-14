# s.napolitanoPhd

This is the repository containing part of the MATLAB code used for PhD Thesis of @SaraNapolitano (from @dibbelab). The other part (regarding the cybergenetics controller to synchronise the cell-cycle) is available at the repository `Cycloop`.


# Repo Contents
+ [TFEB_QuantitativeAnalysis](./TFEB_QuantitativeAnalysis)


# System Requirements
## Hardware Requirements
It requires only a standard computer with enough RAM to run the code in the MATLAB environment.


## Software Dependencies
+ MATLAB R2019a
+ Additional MATLAB toolboxes (e.g. Signal Processing Toolbox, Image Processing Toolbox, etc.)


# Installation Guide
1. Download the repository from https://github.com/dibbelab/s.napolitanoPhd.git .
2. Unpack the files.
3. Start MATLAB and navigate to the `TFEB_QuantitativeAnalysis` folder.


# Instructions for Use
## Fitting
Set the working directory to `./TFEB_QuantitativeAnalysis/Fitting/`. 

+ Run the script `Main_Fitting.m` to derive the values of the time constant used for mTOR inhibition model and reported in Table 4.1. Moreover it generates also the simulations shown in Fig. 4.14B and Fig. 4.14C, i.e. the simulation with the biological hypothesis on the synergistic effect.


## Simulator
Set the working directory to `./TFEB_QuantitativeAnalysis/Simulator/`. 

+ Run the script `Main_Simulator.m` to generate the simulations shown in Fig. 4.6B and Fig. 4.7B, i.e. the simulation without and with the feedback hypothesis, and the simulations shown in Fig 4.15 A-B, i.e. the simulation of TFEB nuclear translocation when the cells are treated with Torin1 in combination with CHX or Bort.


## ExperimentAnalysis
### Segmentation
The folder `./TFEB_QuantitativeAnalysis/ExperimentsAnalysis/Segmentation/` contains the scripts used for segmentation of human HeLa cells. Note that it is needed to run FastER (https://bsse.ethz.ch/csd/software/faster.html) segmentation before the MATLAB script.


### MakeFigure
The folder `./TFEB_QuantitativeAnalysis/ExperimentsAnalysis/MakeFigure/` contains the scripts used to clear the data-set and to generate the images with experimental data shown in the Chapter 4 of the PhD Thesis
