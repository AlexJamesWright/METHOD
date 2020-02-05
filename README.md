<pre><sub><sup>
MMMMMMMM               MMMMMMMMEEEEEEEEEEEEEEEEEEEEEETTTTTTTTTTTTTTTTTTTTTTTHHHHHHHHH     HHHHHHHHH     OOOOOOOOO     DDDDDDDDDDDDD
M:::::::M             M:::::::ME::::::::::::::::::::ET:::::::::::::::::::::TH:::::::H     H:::::::H   OO:::::::::OO   D::::::::::::DDD
M::::::::M           M::::::::ME::::::::::::::::::::ET:::::::::::::::::::::TH:::::::H     H:::::::H OO:::::::::::::OO D:::::::::::::::DD
M:::::::::M         M:::::::::MEE::::::EEEEEEEEE::::ET:::::TT:::::::TT:::::THH::::::H     H::::::HHO:::::::OOO:::::::ODDD:::::DDDDD:::::D
M::::::::::M       M::::::::::M  E:::::E       EEEEEETTTTTT  T:::::T  TTTTTT  H:::::H     H:::::H  O::::::O   O::::::O  D:::::D    D:::::D
M:::::::::::M     M:::::::::::M  E:::::E                     T:::::T          H:::::H     H:::::H  O:::::O     O:::::O  D:::::D     D:::::D
M:::::::M::::M   M::::M:::::::M  E::::::EEEEEEEEEE           T:::::T          H::::::HHHHH::::::H  O:::::O     O:::::O  D:::::D     D:::::D
M::::::M M::::M M::::M M::::::M  E:::::::::::::::E           T:::::T          H:::::::::::::::::H  O:::::O     O:::::O  D:::::D     D:::::D
M::::::M  M::::M::::M  M::::::M  E:::::::::::::::E           T:::::T          H:::::::::::::::::H  O:::::O     O:::::O  D:::::D     D:::::D
M::::::M   M:::::::M   M::::::M  E::::::EEEEEEEEEE           T:::::T          H::::::HHHHH::::::H  O:::::O     O:::::O  D:::::D     D:::::D
M::::::M    M:::::M    M::::::M  E:::::E                     T:::::T          H:::::H     H:::::H  O:::::O     O:::::O  D:::::D     D:::::D
M::::::M     MMMMM     M::::::M  E:::::E       EEEEEE        T:::::T          H:::::H     H:::::H  O::::::O   O::::::O  D:::::D    D:::::D
M::::::M               M::::::MEE::::::EEEEEEEE:::::E      TT:::::::TT      HH::::::H     H::::::HHO:::::::OOO:::::::ODDD:::::DDDDD:::::D
M::::::M               M::::::ME::::::::::::::::::::E      T:::::::::T      H:::::::H     H:::::::H OO:::::::::::::OO D:::::::::::::::DD
M::::::M               M::::::ME::::::::::::::::::::E      T:::::::::T      H:::::::H     H:::::::H   OO:::::::::OO   D::::::::::::DDD
MMMMMMMM               MMMMMMMMEEEEEEEEEEEEEEEEEEEEEE      TTTTTTTTTTT      HHHHHHHHH     HHHHHHHHH     OOOOOOOOO     DDDDDDDDDDDDD <br></sup></sub>
</pre>

[![Build Status](https://travis-ci.org/AlexJamesWright/METHOD.svg?branch=master)](https://travis-ci.org/AlexJamesWright/METHOD)
[![Coverage Status](https://coveralls.io/repos/github/AlexJamesWright/METHOD/badge.svg?branch=master)](https://coveralls.io/github/AlexJamesWright/METHOD?branch=master)
[![Documentation Status](https://readthedocs.org/projects/method/badge/?version=latest)](https://method.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/105871037.svg)](http://doi.org/10.5281/zenodo.1404697)

![alt text](https://github.com/AlexJamesWright/METHOD/blob/master/HighResKHI.gif "High resolution Kelvin-Helmholtz instability")


# Multifluid Electromagneto-HydroDynamics

This is METHOD, a relativistic multi-dimensional, multi-fluid ElectroMagnetoHydroDynamic
solver built and maintained by [Alex Wright](http://cmg.soton.ac.uk/people/ajw1e16/)
under the guidance of [Dr Ian Hawke](https://www.southampton.ac.uk/maths/about/staff/ih3.page).
METHOD is being developed as a result of a PhD project investigating the models of
MHD used to model neutron star mergers. As a result, ideal and resistive single
fluid models exist that are conventional in numerical astrophysics, and a two-fluid
model adapted from Amano 2016.

It also includes the REGIME extension (Wright & Hawke 2019) to ideal MHD. REGIME
allow ideal MHD simulations to capture resistive phenomena by including an
additional source term into the equations of motion.

---------------------------------------------
---------------------------------------------
<br> <br>


## Getting started
To begin using METHOD, first clone the repository

    git clone https://github.com/AlexJamesWright/METHOD.git

To set up the correct directories for storing data, run the provided shell script from the project root,

    bash makePaths.sh

Next, you will need to ensure the Makefiles are valid for your system, changing any compilers to your preferred ones and setting the GoogleTest home directory to its location on you machine. That should be it. Should be.

---------------------------------------------
---------------------------------------------
<br> <br>


## Testing
Once METHOD is installed, check the latest build is working by running the unittests.

We use the Google Test framework for unit testing---any tests are saved in the `Tests/Serial/Src` or `Tests/Parallel/Src` directory. You will need to set the `GTEST_DIR` environment variable (in the Makefile within `Tests/Parallel/Src` and `Tests/Serial/Src`) to point to the GoogleTest root directory.

The serial and parallel versions have separate testing directories. As far as possible the tests are the same, but there is additional testing such that the parallel results match the serial to within floating point accuracy. First, run

    make test

from the `Tests/Serial` directory, then the `Tests/Parallel` directory. To check if results match, from `Tests/Parallel` run

     py.test -v Src/compareSerialAndParallel.py

NOTE: this final test will only pass if the `MATCH_SERIAL` defined constant in `Project/Parallel/Include/timeInt.h` is set to unity, `MATCH_SERIAL=1`.

It is a good idea to check that the examples run successfully next.

---------------------------------------------
---------------------------------------------
<br> <br>


## Example Simulations == Best way to understand functionality!

Example simulations have been provided that illustrate how to use the
various classes. By typing

    make run

in one of the example directories, the relevant object files will be built and
executed. Data is saved in the *Examples/Data* directory and is easily viewed
using the interactivePlot script. We suggest that the spyder environment from
Anaconda is the best way to use this tool.

### InteractivePlot
-------------------
To use the plotting tools, run from the root Example directory something like

    spyder interactivePlot.py

Running this script as main will load any simulation data into the `Plot` object.

This object has a number of pre-made plotting routines, all of which begin with

    Plot.plot

If you wish to create your own plot, you can access the simulation data using the

    Plot.prims

    Plot.cons

arrays, etc. The first index is the variable, followed by `x`, `y`, and `z` index.

To view the available primitive variables, use `Plot.cleanPrimLabels`.

In addition, the Plot object contains a dictionary linking to all the system
constants. For example, to get the value for the adiabatic index used in the
simulation use `Plot.c['gamma']`.

### Animation
-------------

For the Kelvin-Helmholtz simulation, running the `animation.py` script will create an
animatation called `Output.gif` in the root Example directory to view (may take up
to ten minutes to run the simulation and make the animation).

Make sure you clean any timeseries data before running the simulation by running

    bash cleanData.sh

from the root Examples/ directory. The variable being animated can be changed
manually.

---------------------------------------------
---------------------------------------------
<br> <br>


## Builds
To build all the elements of the programme at once go from the Project directory, to either Serial (if you dont have CUDA capable hardware) or Parallel (if you do) and use

    make build

to build each element of the programme.

---------------------------------------------
---------------------------------------------
<br> <br>


## Documentation
I have tried to maintain good documentation standards, but cant guarantee that everything you want will be here. If you are unsure of any functionality, first look in the respective header file, source code, and you are still unsure contact the authors.

The documentation is built automatically using [ReadTheDocs](https://method.readthedocs.io/en/latest/index.html) and can be viewed there.

Alternatively, if you are unsure of any of the functionality, find the documentation in the `index.html` file, generated by [doxygen](https://github.com/doxygen/doxygen) and located in the `Doxumentation/html` directory.
To build the documentation locally simply go the the `Doxumentation` folder and run

    doxygen

## Rootfinder
Some simulations will require the use of an N-dimensional footfinder, either for a (semi-) implicit time integrator or
for the conservative to primitive transformation. We have elected to use the [CMINPACK library](https://github.com/devernay/cminpack)\*, and to use or implement any changes in the library, *cd* into the Cminpack directory and hit

    make objects

to compile all the object files. Then, if the build was successful, don't touch/look at this library again.

---------------------------------------------
---------------------------------------------
<br> <br>

## Simulations

Simulations are run from the *main.cc/cu* scripts. Simply use

    make run

to compile and run the simulation from within `Project/Serial` or `Project/Parallel`. The executable is labelled `main` so

    make build
    ./main

will also work.

---------------------------------------------
---------------------------------------------

<br> <br>

## Saving simulation data

The *Src* directory has a tool for interactively plotting the end state of a simulation. The `interactivePlot.py` script requires data to be saved after the simulation in the *Data*
folder. This is done using the SaveData class---call the class constructor with a pointer to the SimData class whose data you wish to save. Then, simply include

    save.saveAll();

in *main* after the simulation has been evolved. Running the python script as main will load and store the data ready for plotting, and the easiest way to interact with the data is in a python environment such as spyder.

There is also the functionality to save time series data. In order to reduce memory requirements, the user must specify the variables they wish to save (names of the variables should match those given as the labels in the model's header file. To save variables, go into `simulation.cc/cu` and change the three conditional blocks to save the variables you want using

     this->save->saveVar('SomeVar', totalNumberOfUserDefinedVars)

NOTE: The second variable must be included and be the number of variables you wish to save at each output.

---------------------------------------------
---------------------------------------------
<br> <br>


## Authors

[Alex Wright](http://cmg.soton.ac.uk/people/ajw1e16/)  Email: a.j.wright@soton.ac.uk

---------------------------------------------
---------------------------------------------
<br> <br>

\* *due to this cryptic package we have moved bits about and re-ordered various headers and includes. Most of the preprocessor stuff has been deleted (using NVIDIA hardware will result in Cminpack reals defaulting to double precision), some functions have been excluded as they're not needed here, and now for any usage we just include the cminpack.h header file (as opposed to including the CUDA scripts directly).*
