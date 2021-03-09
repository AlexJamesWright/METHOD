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

# Profiling Branch

## Profiling Methods:

* Code profiled is frozen, located at commit: d8c94fd116cd9667147e875ef654f88ed75d63f2, branch: profiling_feb_2021
* Profiled the GPU version using the main script in Project/GPU
* Sim details: 2048x2048, SRMHD, ParallelFlow, KHInstabilitySingleFluid, RK2
* System used: Iridis 5 (GTX1080Ti GPUs, 40 core Intel Skylake CPUs)
* Profiler used: nvvp, nvprof

## Profiling Results:

Using 4 nodes, with 1 MPI process per node, 1 GPU per node and either 1 or 36 OMP threads per MPI process

Memory requirements:
* Each GPU has 12GB memory, each CPU has 192 GB, so need roughly ~10x GPUs to fit the same problem size. For very large problems where the main constraint is fitting in memory rather than time, CPUs likely better than GPUs. 
* 2048x2048 fits on 4 GPUs, 4096x4094 does not

Timestep summary:
* The majority of time in a RK2 timestep is spent in CPU 'useful work', ie not MPI communication
* Neither MPI communication between CPUs, PCI communication between GPU/CPU, nor GPU compute time is significant in relation to CPU compute time. 
* The OpenMP parallelisation of these problem areas on the CPU seems to be inefficient -- in the timestep breakdown below it looks like going from 1 thread to 36 threads speeds the cpu code by only roughly 3x. This is bourne out by the timing of the non-GPU code, which is fully MPI parallelised and is quicker than the GPU version. It may be that the CPU node is being shared between multiple GPU jobs for different users, or the OMP parallelisation may be able to be improved. 

* Conclusions: 
1) further GPU optimisation work should focus on moving more parts of the timestep to the GPU, without worrying too much about the cost of copying data between GPU/CPU at least at first
2) enabling OpenMP for the GPU version is important. This will parallelise the problem areas on the CPU and can be enabled from current makefiles. Will also need to `export OMP_NUM_THREADS=[threads]` in job submission script. 
3) It worth looking into and improving the efficiency of the OMP parallelisation if the work cannot be instead moved to the GPU

Timestep breakdown:

Run A: 4 nodes, each running 1 proc x 1 OMP thread (ie most CPU cores are unused). Timescale: 1 RK2 timestep = 10.31 seconds

Total time for 21 timesteps: 224.44s

See [RK2.cu](Project/GPU/Src/RK2.cu) for exact details on how exactly the code is grouped and annotated to show the different components of the timestep in alternating green and black. In order, the sections are: cons2prims, flux method, cons2prims 2, boundary condition calculation, flux method 2, construct solution, get prims, boundary condition calculation 2.

![rk2_omp1](rk2_omp1.png)

Run B: 4 nodes, each running 1 proc x 36 OMP threads. Timescale: 1 RK2 timestep = 3.02 seconds

Total time for 21 timesteps: 65.19s

![rk2_omp36](rk2_omp36.png)

Run C: 4 nodes, each running 36 MPI procs x 1 OMP threads, no GPU. 

Total time for 21 timesteps: 8.14s


---------------------------------------------
---------------------------------------------
<br> <br>

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
