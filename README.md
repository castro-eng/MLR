# MLR

### What is this repository for?

MLR is a software developed for simulation of reservoirs containing multilayer laminar structures, such as faults and fractures.

### How do I get set up?

Follow the instructions:

1. Copy the repository folder (*** Downloads *** from the menu on the left of this page) in a Linux distribution.

2. Have gfortran and h5fc installed in the machine.

### Usage

There are two scripts for the compilation and usage of the program.

The file "compilar.sh" is a script for compiling the fortran code.

The file "rodarSimulado.sh" is a script for running the MRL program, and it's necessary to indicate the path of the experiment to be simulated. (Ex-> ./rodarSimulador.sh experiment/Exp1/HighFidelity01 )

OUTPUT FILES (in the "out" folder of the experiment):

- tela.txt --> A summary of the main simulation information. Containing at the end of the file the equivalent permeability and transmissibility information.
- ExperimentoBloco.xdmf --> a paraview file format in XDMF with the results for the simulated model.
- <name>\_pic.vtk --> several paraview file from the simulated model.

### Who do I talk to?

Eduardo da Silva Castro, National Laboratory of Scientific Computing
**_ castro.eng@gmail.com _**

### License

Copyright (C) 2022 Eduardo da Silva Castro
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
