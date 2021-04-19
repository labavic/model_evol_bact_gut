This code contains a Fortran 90 implementation of the code described in "Modeling evolution of bacteria in the gut", by Darka Labavic, Claude Loverdo, and Anne-Florence Bitbol.

Archived version:

We solve a set of coupled partial differential equations that describe the dynamics of food, wild type and mutant bacteria concentrations. The code first solves PDEs without mutant bacteria in the system for a given fixed time. Second, it takes the final state from step one and solves the system with an initial amount of mutants, considering several possible mutant initial positions.

The output gives three text files:

parameters.dat - contains all parameters used in the run

Final_Stat.dat - stationary state of food and wild type (result of step one) 

dataMutantEnd.dat - state of the system after the second step. Each row corresponds to one initial position. There are 3N+1 columns, the first one corresponds to the initial mutant position, the next N columns to wild type concentrations at each of the N discrete spatial points, then the following N to mutant concentrations, and final N columns to food concentrations.
The code can be compiled by executing the command: gfortran MBEG.f90, and then executed by: ./a.out

Parameters can be changed directly in the MBEG.f90 file (note that the code then needs to be recompiled) or they can be given in the command line as follows: ./a.out  par_chr1 par_val1 par_chr2 par_val2 etc. 

where par_chr is a parameter name (listed below) and par_val is the corresponding parameter value.

d - diffusion coefficient

v - flow velocity

a - conversion parameter between food and bacteria units

k - Monod constant

r - growth rate of the wild type

rm - growth rate of the mutant 

Fin - food inflow concentration

dx - discretization step in space

dt - discretization step in time

time - total time to integrate in step one and in step two

len - length of the space segment

mi - 0. if the initial amount of mutant is the same for all initial positions, 1. if it is proportional to the reproduction rate at the initial position

M0 - initial local amount of mutant

The variable "path" in the code needs to be changed to the desired path for the output.

The source code is freely available under the GNU GPLv3 license.
