This code contains a Fortran 90 implementation of the code described in "Modelling evolution of bacteria in the gut", by darka Labavic, Claude Loverdo, and  Anne-Florence Bitbol.
Archived version:

We solve a set of coupled partial differential equations that describe dynamics of the interaction of food, wild type and mutant bacteria. The code first solves PDEs without mutant bacteria in the system for a fixed given time. Second, it takes the state from step two and solves the system with an initial amount of the mutant at several positions in space, one at the time.

The output gives three text files:

parameters.dat - contains all parameters used in the run

Final SS.dat  - stationary state of food and wild type (result of step one)

dataMutantEnd.dat - state of the system after second step. Each row corresponds to the one initial position. There are 3N+1 columns, first corresponds to the initial position, nex N columns to the state of the wild type, second N to the state of mutant, and third N to the state of food.

The code can be compiled by command: gfortran MBEG.f90, and then executed by: ./a.out

Parameters can be changed directly in the MBEG.f90 file (code needs to be recompiled if changed like that) or they can be given in the command line after like ./a.out  par_chr par_val 

where par_chr is a parameter name (listed below) and par_val parameter value.

d - diffusion coefficient

v - flow velocity

a - conversion parameter between food and bacteria units

k - monod constant

r - growth rate of the wild type

rm - growth rate of the mutant

Fin - food inflow concentration

dx - discretization step in space

dt - discretization step in time

time - total time to integrate in step one and in step two 

len - length of the space segment

mi - 0. if the initial amount of mutant is the same for all initial positions, 1. if it is proportional to the reproduction rate at the initial position

M0 - initial amount of the mutant 

The variable "path" in the code needs to be changed to the desired path for the output.

The source code is freely available under the GNU GPLv3 license.


