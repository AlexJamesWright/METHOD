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

![alt text](https://github.com/AlexJamesWright/METHOD/blob/master/METHODAdvert.gif "METHOD Advert: kudos to pyro2 for the inspo")


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

### Quick start
To begin using METHOD, first clone the repository

    git clone https://github.com/AlexJamesWright/METHOD.git

To set up the correct directories for storing data, run the provided shell script from the project root,

    bash makePaths.sh

Compile a simple example: single fluid, Kelvin-Helmholtz instability with random perturbation. This example will run in serial on one CPU core and should only require the gnu c++ compiler g++ to be installed. 

```
cd Examples/KelvinHelmholtz/SingleFluidIdealSerial/
make
```

Run the example using:
```
./main 
```

This will run a small number of timesteps and save the final simulation state in plaintext form in Examples/Data/Final. 

```
t = 0.000000
t = 0.167385
t = 0.498842
t = 0.830298
t = 1.161754
t = 1.493211
t = 1.824667
t = 2.156123
t = 2.487580
t = 2.819036
```

For instructions on running other simulations on a range of computer architectures and with different options for data output format, see [workflows.md](workflows.md).

### Running tests

The following instructions will run a subset of unit tests for testing the core functionality of METHOD when running in serial on a single CPU core and outputting data in plain text format. These tests require the gnu c++ compiler g++ and Python 3 to be installed. For instructions on running the full test suite, which tests a range of computer architectures and output formats, see [workflows.md](workflows.md).  

First, make sure to have cloned METHOD and run `bash makePaths.sh` as above. 

Set up a python virtual environment with modules required for testing by typing the following in the root directory:

```
python3 -m venv venv
source venv/bin/activate
python -m pip install -r Scripts/IridisEnv/requirements.txt
```

Clone the GoogleTest repository into the directory above the METHOD root directory:

```
git clone https://github.com/google/googletest.git
```

Run the serial CPU plain text tests using:

```
cd Tests/CPU
make test_serial
```

This will run a number of tests, and should end with all tests passing in output similar to:

```
...
[----------] Global test environment tear-down
[==========] 21 tests from 5 test suites ran. (7238 ms total)
[  PASSED  ] 21 tests.
```

---------------------------------------------
---------------------------------------------
<br> <br>

## Documentation
I have tried to maintain good documentation standards, but cant guarantee that everything you want will be here. If you are unsure of any functionality, first look in the respective header file, source code, and you are still unsure contact the authors.

The documentation is built automatically using [ReadTheDocs](https://method.readthedocs.io/en/latest/index.html) and can be viewed there.

Alternatively, if you are unsure of any of the functionality, find the documentation in the `index.html` file, generated by [doxygen](https://github.com/doxygen/doxygen) and located in the `Doxumentation/html` directory.
To build the documentation locally simply go the the `Doxumentation` folder and run

    doxygen
    
---------------------------------------------
---------------------------------------------
<br> <br>

## Authors

[Alex Wright](http://cmg.soton.ac.uk/people/ajw1e16/)  Email: a.j.wright@soton.ac.uk <br>
[Ania Brown](https://github.com/aniabrown) <br>
[Sam Mangham](https://github.com/smangham) <br>
[Ian Hawke](https://cmg.soton.ac.uk/people/ih3/)
---------------------------------------------
---------------------------------------------
<br> <br>
