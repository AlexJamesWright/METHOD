
@mainpage
@section intro Multi-Fluid ElectroMagnetoHydroDynamics
@par
This is METHOD, a relativistic multi-dimensional, multi-fluid ElectroMagnetoHydroDynamic
solver built and maintained by [Alex Wright](http://cmg.soton.ac.uk/people/ajw1e16/)
under the guidance of [Dr Ian Hawke](https://www.southampton.ac.uk/maths/about/staff/ih3.page).
METHOD is being developed as a result of a PhD project to model neutron star systems
with multi-fluid models. As a result, ideal and resistive single fluid models exist
that are more conventional in astrophysical model, and a relatively new two-fluid
model adapted from Amano 2016.

@section install Installation
@par
The easiest way to install METHOD is to clone the github [repository](https://github.com/AlexJamesWright/METHOD).
The repository contains all the source files, tests and makefiles that one needs
to get going. To start, copy the HTTPS tag on the main page or run

    git clone https://github.com/AlexJamesWright/METHOD
to copy all files to your machine. There are then multiple options, but the safest
is to go to the *Project* directory and build all objects first. Do this by running

    make build
This will compile all the source code and build the rootfinders, and also generate
the documetation. To check the build has worked, the best option is to run all the
tests. From the *Test* directory, running

    make test
will compile all unit tests and, using the GoogleTest framework, execute the
available tests. If any tests fail unexpectedly, update the files on you machine
by pulling the repository and try again.

@section example Example Simulations
@par
There is currently work on going to build a set of test cases and examples from  
which an understanding of the underlying code and processes that are implemented
can be understood. These examples range from standard Brio Wu shock tube tests
to set ups in which there are exact solutions in the form of self similar current
sheets. These examples should be the first port of call for anyone trying to use
the METHOD code for the first time, of if you are using a particular aspect of the
code for the first time and wish to understand what is required to perform a succesfull
simulation.
