# BrownianBead
*Author : Martin Rieu, http://www.normalesup.org/~mrieu*
*March 2021* 

## Description 

C++ code for the simulation of a brownian magnetic bead tethered to a dsDNA close to a surface. Can be used to simulate the expected noise for superresolution magnetic tweezers experiments (see for example _Rieu, Martin, et al. "Parallel, linear, and subnanometric 3D tracking of microparticles with Stereo Darkfield Interferometry." Science Advances 7.6 (2021): eabe3902._)

The code integrates the hydrodynamic corrections due to the presence of a wall as summarized in _Series for computation of transverse translationnal and rotationnal coefficients. Physica A 189 (1992) 447-477  G.S. Perkins and R.B. Jones_). 

The code can be used as a plug-in for the software _Xvin_ that can be found in the repository https://tig.phys.ens.fr/ABCDLab/xvin/. Once the plug-in is loaded, it can be called through the menu "Treatment -> BrownianBeads -> Simulate Brownian Bead"

It can also be used as a stand-alone software but in this case, no graphical interface is provided. 
The software can be easily compiled

## Dependencies
-Eigen3
-gsl

## Linux installation of the dependencies
sudo apt install libeigen3-dev libgsl-dev

## Linux command for compilation 
g++  HydroCorrections.cpp BrSimulator.cpp  -o BrownianBead \`pkg-config --cflags --libs eigen3\` -lgsl -O3

## Stand-alone output
The stand-alone function has been written for debug purposes. 
If you wish to extract data from the simulation, you need to modify the function main() in the file BrSimulator.cpp and write the data stored in the float arrayx named x,y, and z in the container of your choice. 
