# MOLTO-IT
**MOLTO-IT** (Multi-Objective Low-Thrust Optimizer for Interplanetary Trajectories) is a fully automated Matlab tool for the preliminary design of low-thrust, multi-gravity assist trajectories. The software combines an outer loop that provides multi-objective optimization via a genetic algorithm (NSGA-II) with an inner loop that supplies gradient-based optimization (fmincon) of a shape-based low-thrust trajectory parameterization.

It includes simplifications such as coplanar bodies and no enforced propulsion constraints along with the shape-based parameterization of the low-thrust arc. As a result, run-times are reduced to the order of minutes, which makes the tool able to run on a standard desktop computer. 

Users only need to input a series of parameters, such as the spacecraft’s departure body, its final destination and some physical characteristics, as well as the range of launch dates, flight times and a list of available planets to flyby. The software tool then uses these data points to automatically compute the set of low-thrust trajectories, including the number, sequence and configuration of flybys that accomplish the mission with the minimum cost. Candidate trajectories are evaluated and compared in terms of total flight time and propellant mass consumed. 

### Goal
The purpose of MOLTO-IT is to provide a fast and robust mission design environment that allows the user to quickly and inexpensively perform trade studies of various mission configurations and conduct low-fidelity analysis.

### Mathematical formulation
For a detailed mathematical description of the methodology, have a look at: 

David Morante, Manuel Sanjurjo Rivo, and Manuel Soler.  "Multi-Objective Low-Thrust Interplanetary Trajectory Optimization Based on Generalized Logarithmic Spirals", Journal of Guidance, Control, and Dynamics, Vol. 42, No. 3 (2019), pp. 476-490. 
https://doi.org/10.2514/1.G003702

## Installation Guide

### Dependencies
fmincon
Spice 
## Quick Usage Guide


