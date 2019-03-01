# MOLTO-IT
[![DOI](https://zenodo.org/badge/169291732.svg)](https://zenodo.org/badge/latestdoi/169291732)

**MOLTO-IT** (Multi-Objective Low-Thrust Optimizer for Interplanetary Trajectories) is a fully automated Matlab tool for the preliminary design of low-thrust, multi-gravity assist trajectories. The software combines an outer loop that provides multi-objective optimization via a genetic algorithm (NSGA-II) with an inner loop that supplies gradient-based optimization (fmincon) of a shape-based low-thrust trajectory parameterization. It includes simplifications such as coplanar bodies and no enforced propulsion constraints along with the shape-based parameterization of the low-thrust arc.

MOLTO-IT is designed to work with minimal user oversight. Users only need to input a series of parameters, such as the spacecraft’s departure body, its final destination and some hardware characteristics (Launcher vehicle, mass, propulsion), as well as the range of launch dates, flight times and a list of available planets to flyby. The software tool then uses these data points to automatically compute the set of low-thrust trajectories, including the number, sequence and configuration of flybys that accomplish the mission most efficiently. Candidate trajectories are evaluated and compared in terms of total flight time and propellant mass consumed. 

MOLTO-IT is easy to use and capable of designing a wide variety of missions on standard desktop computer.Computational times ranges from the order of minutes to a couple of hours, depending on the size of ther search spacee.

### Goal
The purpose of MOLTO-IT is to provide a fast and robust mission design environment that allows the user to quickly and inexpensively perform trade studies of various mission configurations and conduct low-fidelity analysis.

### Mathematical formulation
For a detailed mathematical description of the methodology, have a look at: 

David Morante, Manuel Sanjurjo Rivo, and Manuel Soler.  "Multi-Objective Low-Thrust Interplanetary Trajectory Optimization Based on Generalized Logarithmic Spirals", Journal of Guidance, Control, and Dynamics, Vol. 42, No. 3 (2019), pp. 476-490. 
https://doi.org/10.2514/1.G003702

## Installation Guide
Installation requires simply that you download [MOLTO-IT](https://github.com/uc3m-aerospace/MOLTO-IT/) and add the base directory to your Matlab path.

### Dependencies
A recent version of Matlab is needed to run the code (R2016a or newer). The [Matlab Optimization Toolbox](https://es.mathworks.com/help/optim/index.html) is required to use most of the functionality.

Download the Matlab interface of the [Spice Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html) for planets, comets and asteroids ephemerides calculations. Include the downloaded folders and subfolders in **MOLTO-IT/spice**. 

Download apropriate ephemerides [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/) .bsp kernels. Include them in **MOLTO-IT/spice/kernels** 

Modify the function **MOLTO-IT/spice/load_spice_kernels** to load only your rewuired kernels.

## Quick Usage Guide


