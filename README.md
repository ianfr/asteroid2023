# Asteroid Simulation
## About
3D rubble-pile asteroid simulation. Tested with Intel Fortran on Ubuntu 22.04 64-bit.

## Goals
- [x] Generate asteroids made of point particles according to the user specification
- [x] Generate data for easy post-processing with Paraview
- [x] Simulate the motion of the particles under gravity
- [x] Perform collision detection
- [ ] Implement cohesion
- [ ] Parallelize for CPU

## Simulation Details:
Forces:
- Gravity

Integration options:
- Euler (order 1)
- Verlet (order 2, **symplectic**) [`WIP`]
- Runge-Kutta (order 4)

Collision schemes:
- Complete brute force (_O(n^2)_ time complexity)
- Radially constrained brute force with probablistic complete brute force (_O(n*R^3 + eps*n^2)_) [`WIP`]

Adaptive time-stepping based on:
- Maximum particle velocity
- Mean particle velocity [`WIP`] 

## Building & Running
This project uses the CMake build system. To build and run, execute the following commands:
* `cd asteroid2022`
* `mkdir build && cd build`
* `cmake ..`
* `make`
* `./prog` to run the program.

## Configuration
The first step for configuring a run is to set high-level parameters in the main file *main_prog.f90* (like the timestep size & number of asteroids).

Additionally, the output directory set in *control.txt* must exist before the program runs. Assuming you're already in the project directory, the command would be `mkdir OUT/X` where `X` is the value of the output directory field in *control.txt*.

## Viewing the results
A Python script is included to facilitate viewing results when Paraview cannot load too many files at one time, but **Paraview is the preferred method**.

Paraview allows for opening "grouped" CSV files in numerical filename order, which is what the program outputs. In Paraview:
* Open the CSV group 'ast.*' in the relevant output directory
* Click 'Apply' and then X out of the default Table View
* Right-click on the group then do Add Filter -> Alphabetical -> Table to Points
* Set the X, Y, and Z Columns to rx, ry, and rx respectively with the dropdown menus
* Hit 'Apply'
* Change 'Representation' to 'Point Gaussian'
* Change 'Coloring' to 'col_ast'

To view a frame-by-frame animation, just click the play button at the top!

**Repository copyright Ian Friedrichs 2023. All rights reserved.**