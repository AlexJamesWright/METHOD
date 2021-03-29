
# METHOD workflows

## Introduction

This document outlines detailed instructions for compiling and running METHOD on a variety of computer architectures. The document steps through installing dependencies, running the example scripts in `Examples`, running the test suite in `Tests`, and finally creating new simulation scripts in the `Project` directory. See the [README](README.md) for instructions on getting started with just the simplest available example. 

METHOD can run serially on a single CPU core and optionally also use OpenMP (suitable for using multiple CPU cores within a single computer) and MPI (suitable for running across several nodes in a compute cluster environment). There is also a GPU version which implements a more limited number of models and initial conditions. These options require a compiler which implements OpenMP 3.1, the MPI libraries (OpenMPI 2.0 or later, or mpich), and an NVIDIA GPU with CUDA drivers and toolkit, respectively. 

Both the CPU and GPU versions of METHOD can output data as either plain text or HDF5 format. To output HDF5 from a single process, either the serial or parallel HDF5 libraries must be installed. To output HDF5 from a multiprocess MPI code, the parallel version of HDF5 must be installed. 

Four target systems will be used as examples in these instructions. 

1. Ubuntu 18.04 (Bionic) with GTX 1050 NVIDIA GPU
2. Ubuntu 20.04 (Focal) with GTX 1050 NVIDIA GPU
3. MacOS 10.15 (Catalina) with no GPU 
4. Iridis 5 supercomputer compute node, RHEL 7.0 with Module system for software, with GTX 1080 NVIDIA GPU

---------------------------------------------
---------------------------------------------
<br> <br>

## Installing METHOD and dependencies 

To begin using METHOD, first clone the repository

    git clone https://github.com/AlexJamesWright/METHOD.git

Then set up the correct directories for storing data by running the provided shell script from the project root,

    bash makePaths.sh

### OpenMP

**Linux**: the GNU (gcc/g++) version 4.8 or later supports OpenMP 3.1 for c/c++. This or later versions are the default on most Linux systems. 

**MacOS Catalina**: By default, gcc maps to Apple Clang 11.0, which does not support OpenMP. gcc can in theory be installed using homebrew, though in practice this I would recommend disabling OpenMP (set USE_OMP=0 in and Makefiles used) for MacOS as MPI gives similar performance and seems to be easier to install.

### MPI

> **Note**: While it is possible to use multiple versions of MPI on the same system, it is easier if there is only one. Check if MPI is already installed using `mpirun --version`. Currently supported versions are MPICH 3.0 or later and OpenMPI 4.0 or later.  

**Ubuntu 20.04**: Install the MPI libraries using `sudo apt-get install libopenmpi-dev` or `sudo apt-get install libmpich-dev` This will install OpenMPI v4 or MPICH v3. 

**Ubuntu 18.04**: Install the MPI libraries using `sudo apt-get install libmpich-dev`. This will install MPICH v3. Note the version of OpenMPI available through the package manager (v2) is not supported by METHOD.

**Iridis 5**: Load the mpich module using `module load mpich/3.2.1/gcc`. Note that this will need to be done in the login node if compiling on the login node, and also in the batch script used to run METHOD. 

**MacOS Catalina**: The two options are to install using [homebrew](https://brew.sh/
) or to build from source. 

1)
Install homebrew using instructions here: https://brew.sh/. Then:
```
brew update
brew install open-mpi
```

2)

```
curl https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.0.tar.gz --output openmpi-4.1.0.tar.gz
gunzip -c openmpi-4.1.0.tar.gz | tar xf -
cd openmpi-4.1.0
./configure --prefix=/usr/local
<...lots of output...>
make all install
```

> **Known issues**: 
OpenMPI v2 has been found to lead to stack smashing errors when creating an HDF5 file in parallel. See Issue [#31](https://github.com/AlexJamesWright/METHOD/issues/31)

### CUDA toolkit and drivers

**Ubuntu**: To install the latest version of CUDA available through the package manager, use `sudo apt install nvidia-cuda-toolkit`. Otherwise, to install a specific version, follow the instructions here: https://developer.nvidia.com/Cuda-Toolkit-archive. METHOD has been tested with CUDA 8 and CUDA 10.

**Iridis 5**: Load the CUDA module using `module load cuda/8.0`. Note that this will need to be done in the login node if compiling on the login node, and also in the batch script used to run METHOD. 

### HDF5

**Ubuntu**: Depending on which version of MPI is installed, install HDF5 using either `sudo apt-get install libhdf5-openmpi` or `sudo apt-get install libhdf5-mpich`

**Iridis 5**: If not using mpi, load the serial hdf5 module using `module load hdf5/1.10.2/gcc/serial`. If using mpi, use `module load hdf5/1.10.2/gcc/parallel`.  Note that this will need to be done in the login node if compiling on the login node, and also in the batch script used to run METHOD. 

**MacOS Catalina**: Install homebrew using instructions here: https://brew.sh/. Then use `brew install hdf5` for the serial version or `brew install hdf5-mpi` for the mpi version.

---------------------------------------------
---------------------------------------------
<br> <br>

## Running example simulations == best way to understand functionality!

Example simulations have been provided that illustrate how to use the various classes. 

> **Note**: Currently only the examples in Examples/KelvinHelmholtz are up to date with the most recent version of METHOD. All examples are configured to run on CPU, not GPU.

There are two initial conditions modelled in this folder: SingleFluidIdeal and SingleFluidIdealRandom. For these initial conditions, the following configurations are available:

* Serial: Run on a single CPU core, output data to plain text format
* SerialHDF5: Run on a single CPU core, output data to hdf5 format
* MPI: Run on multiple CPU cores using MPI, output data to plain text format
* MPIHDF5: Run on multiple CPU cores using MPI, output data to plain text format

There are also SingleFluidIdealRandomMPIHDF5CheckpointRestart and SingleFluidIdealRandomSerialHDF5CheckpointRestart, which illustrate checkpoint restart by initialising the simulation from a restart file written out during a previous SingleFluidIdealRandom simulation (in either serial or parallel).

### Building examples

In general, you should be able to build the examples simply by typing `make` in each example directory (make sure to have loaded the correct modules if using Iridis 5). 

There are some instances where you will need to edit the first part of the Makefile, shown below. If your system has multiple versions of MPI running, you may need to specify the particular version of the mpi compiler to use through the MPI_CC variable. If using HDF5, also make sure that the correct version of MPI is being used by the HDF5 compiler using the HDF5_CC variable. 

```
# -------------- PARAMETERS USERS ARE LIKELY TO NEED TO EDIT -------------------
# Whether to use MPI for multi-cpu processing
USE_MPI = 1
USE_OMP = 0
USE_HDF = 1

# Compiler
CC = g++
# --- if USE_MPI ---
# If using mpi, aditionally specify a c++ capable mpi compiler. In systems with multiple versions of MPI, 
# the particular version may need to be specified with eg mpicxx.mpich
MPI_CC = mpic++
# --- if USE_HDF ---
# If using hdf5, additionally specify a hdf5 compiler. If using mpi, this must be the version of the hdf5 
# compiler available on your system that links the correct mpi libraries. Should 
# be one of h5pcc, h5pcc.openmpi or h5pcc.mpich. 
HDF5_CC = h5pcc

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------
```

### Running examples locally after building

**Serial**: To run a serial examples locally, type either `./main` or `./main [seed]` to run the SingleFluidIdeal or SingleFluidIdealRandom examples respectively.

**Parallel**: To run examples which use MPI to launch multiple processes, type `mpirun -np 4 ./main [seed]`. This will run the code using 4 processes, the process count expected by the example codes. 

### Running examples in a batch environment (eg Iridis 5)

To run any of the MPI examples using Iridis 5, submit the example job located in the Scripts/IridisEnv directory from the example directory:

```
cd Examples/KelvinHelmholtz/SingleFluidIdealRandomMPIHDF5
sbatch ../../../Scripts/IridisEnv/examples_cpu_job.sh
```

This will load the correct modules, do a clean build of the example and run it in parallel. 

You can also test building the example from the login node by load the modules from `example_job.sh` and typing `make`. 

### Viewing data

Data is saved in the *Examples/Data* directory 

> **Note**: The InteractivePlot and Animation tools are no longer supported: Use at your own risk.

#### InteractivePlot

The example data can be viewed using the interactivePlot script. We suggest that the spyder environment from
Anaconda is the best way to use this tool.

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

#### Animation

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


## Running tests

Once METHOD is installed, check the latest build is working by running the unit tests. 

We use the Google Test framework for unit testing---any tests are saved in the `Tests/CPU/Src` or `Tests/GPU/Src` directory. 

As far as possible the tests in the two testing directories are the same, though the full simulation tests are less comprehensive on GPU than on CPU. For both the CPU and GPU code, the core unit tests do not use MPI, but there is additional testing of full simulation runs such that the serial and MPI results match to within floating point accuracy. The outputs of CPU and GPU runs are not compared as the different order of floating point calculations between the CPU and GPU versions causes the results to diverge for some simulations. 

The CPU test suite additionally contains tests that compare the outputs of full simulations in the older plain text format with outputs in the newer HDF5 format, to verify the correctness of the HDF5 writer.  

To run tests, first set up a python virtual environment with modules required for testing by typing the following in the root directory:

```
python3 -m venv venv
source venv/bin/activate
python -m pip install -r Scripts/IridisEnv/requirements.txt
```

Clone the GoogleTest repository into the directory above the METHOD root directory:

```
git clone https://github.com/google/googletest.git
```

### CPU tests

To run the full CPU test suite, the MPI and HDF5 dependencies need to be installed. Then, run

```
cd Tests/CPU
make test
```

On Iridis, instead submit this work as a job from the Tests/CPU folder using `sbatch ../../Scripts/IridisEnv/tests_cpu_job.sh`

This will run a number of tests. Note that these tests are currently in separate batches that are not all summarised at the end of the test output, so it is currently necessary to look through the whole job output for any failing tests. 

If neither MPI nor HDF5 is installed, a smaller core of unit tests can be run using:

```
cd Tests/CPU
make test_serial
```

### GPU tests

The run the full GPU test suite, the MPI, HDF5 and CUDA dependencies need to be installed. Then, run

```
cd Tests/GPU
make test
```

On Iridis, instead submit this work as a job from the Tests/GPU folder using `sbatch ../../Scripts/IridisEnv/tests_gpu_job.sh`

As for the CPU version, these tests are currently in separate batches that are not all summarised at the end of the test output, so it is currently necessary to look through the whole job output for any failing tests. 

---------------------------------------------
---------------------------------------------
<br> <br>



## Creating your own simulation scripts

To create your own METHOD script, edit main.cc in Project/CPU/Src or main.cu in Project/GPU/Src. 

### Compiling and running CPU scripts

Edit the first part of the Makefile to reflect whether you wish to use MPI, OpenMP and/or HDF5, and to specify your mpi and hdf5 compilers:

```
# -------------- PARAMETERS FOR USERS TO EDIT --------------------

# Whether to use MPI for multi-cpu processing
USE_MPI = 1
USE_OMP = 0
USE_HDF = 1

# Compiler
CC = g++
# --- if USE_MPI ---
# If using mpi, aditionally specify a c++ capable mpi compiler. In systems with multiple versions of MPI, 
# the particular version may need to be specified with eg mpicxx.mpich
MPI_CC = mpic++
# --- if USE_HDF ---
# If using hdf5, additionally specify a hdf5 compiler. If using mpi, this must be the version of the hdf5 
# compiler available on your system that links the correct mpi libraries. Should 
# be one of h5pcc, h5pcc.openmpi or h5pcc.mpich. 
HDF5_CC = h5pcc

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------
```

Then compile with

    make 

**Serial**: Run from Project/CPU using `./main` 

**Parallel**: Run using `mpirun -np [nprocs] ./main`. Where nprocs is the number of processes to use. Note that nprocs must equal nxRanks x nyRanks x nzRanks in `ParallelEnv env(argc, argv, nxRanks, nyRanks, nzRanks)`. To run on Iridis 5, edit `Scripts/IridisEnv/examples_cpu_job.sh` so that --ntasks-per-node x nodes equals nprocs and replace `mpirun -np 4` with `mpirun -np [nprocs]`, then run from Project/CPU using `sbatch ../../Scripts/IridisEnv/examples_cpu_job.sh`. 

### Compiling and running GPU scripts

Edit the first part of the Makefile to reflect whether you wish to use MPI, and/or HDF5, and to specify your mpi compiler and path to hdf5 libraries. Unfortunately we can't simply use the hdf5 compiler wrapper here as it seems to interact badly with the gpu compiler nvcc. 

Importantlly, make sure to set GPU_COMPUTE_CAPABILITY correctly according to the list here: https://developer.nvidia.com/cuda-GPUs (remove the decimal point, so cc 5.2 would be 52). You can find your GPU model by typing `nvidia-smi`. If compute capability is wrong, the code will likely still compile but may lead to wrong answers. 

```
# -------------- PARAMETERS FOR USERS TO EDIT --------------------

# if USE_MPI=1, need to use parallel versions of objects, such as ParallelEnv, ParallelSaveData etc
USE_MPI=1
USE_HDF=1

# The compute capability of the GPU 
GPU_COMPUTE_CAPABILITY = 52

# --- IF USE_MPI ---
# The c++ capable mpi compiler. In systems with multiple versions of MPI, the particular version may need to be specified with eg
# mpicxx.mpich
MPI_CC = mpic++

# --- IF USE_HDF ---
# HDF5 libraries must be linked explicitly like this rather than using the hdf5 compiler h5pcc. 
# h5pcc should wrap mpicc with the hdf5 libraries included, but appears to interact badly with nvcc
# The library paths below are found using h5pcc -show
HDF5_FLAGS = -I/local/software/szip/2.1.1/include -L/local/software/hdf5/1.10.2/gcc/parallel/lib -L/local/software/szip/2.1.1/lib -lsz -lz -ldl -lm -I/local/software/hdf5/1.10.2/gcc/parallel/include -lhdf5 -lhdf5_hl
# Ubuntu 18.04 mpich example
#HDF5_FLAGS = -I/usr/include/hdf5/mpich -L/usr/lib/x86_64-linux-gnu/hdf5/mpich /usr/lib/x86_64-linux-gnu/hdf5/mpich/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/mpich/libhdf5.a -lsz -lz -lm

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------
```

Then compile with

    make 

**Serial**: Run from Project/CPU using `./main` 

**Parallel**: Run using `mpirun -np [nprocs] ./main`. Where nprocs is the number of processes to use. Note that nprocs must equal nxRanks x nyRanks x nzRanks in `ParallelEnv env(argc, argv, nxRanks, nyRanks, nzRanks)`. To run on Iridis 5, edit `Scripts/IridisEnv/examples_gpu_job.sh` so that --ntasks-per-node x nodes equals nprocs and replace `mpirun -np 4` with `mpirun -np [nprocs]`, then run from Project/CPU using `sbatch ../../Scripts/IridisEnv/examples_gpu_job.sh`. 


### Saving simulation data

The *Src* directory has a tool for interactively plotting the end state of a simulation. The `interactivePlot.py` script requires data to be saved after the simulation in the *Data*
folder. This is done using the SaveData class---call the class constructor with a pointer to the SimData class whose data you wish to save. Then, simply include

    save.saveAll();

in *main* after the simulation has been evolved. Running the python script as main will load and store the data ready for plotting, and the easiest way to interact with the data is in a python environment such as spyder.

There is also the functionality to save time series data. In order to reduce memory requirements, the user must specify the variables they wish to save (names of the variables should match those given as the labels in the model's header file. To save variables, go into `simulation.cc/cu` and change the three conditional blocks to save the variables you want using

     this->save->saveVar('SomeVar', totalNumberOfUserDefinedVars)

NOTE: The second variable must be included and be the number of variables you wish to save at each output.


## Rootfinder
Some simulations will require the use of an N-dimensional footfinder, either for a (semi-) implicit time integrator or
for the conservative to primitive transformation. We have elected to use the [CMINPACK library](https://github.com/devernay/cminpack)\*. This library is built automatically from the build scripts used in Project, Examples and Tests and should not need to be touched. However, if you need to implement and test any changes to the library, you can also *cd* into the Cminpack directory and hit

    make objects

to compile all the object files. 

---------------------------------------------

\* *due to this cryptic package we have moved bits about and re-ordered various headers and includes. Most of the preprocessor stuff has been deleted (using NVIDIA hardware will result in Cminpack reals defaulting to double precision), some functions have been excluded as they're not needed here, and now for any usage we just include the cminpack.h header file (as opposed to including the CUDA scripts directly).*

---------------------------------------------
---------------------------------------------
<br> <br>




