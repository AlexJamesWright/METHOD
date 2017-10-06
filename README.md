# MultiDMHD
-----------------------
------------------------

**M**ulti-**D**imensional **M**agneto**H**ydro**D**ynamics
---------------------------------------------

The extension of the 2DEuler python code, based in CUDA and C++.

## Testing
We use the Google Test framework for unit testing. *Dont touch the GoogleTest directory!* Any tests are saved in the `Tests/Src` directory.

To build and run tests, from the `Tests` directory and run

  `make tests`
  
Any failed tests will be highlighted in red, and if you bbroke it, you gotta fix it.

## Simulations
Simulations are run from the *main.cu* scripts. The Makefile in the Project folder links to the test directory, it is recommended that simulations are run with 

  `make all`
  
so that gtest will flag up any broken tests since the last change. Otherwise, simply use
  
  `make run`
    
to compile and run the simulation.

Happy modelling!
