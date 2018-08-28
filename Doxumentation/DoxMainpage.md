
@mainpage
@section intro Multi-Fluid ElectroMagnetoHydroDynamics
@par
<pre><br><br>
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
MMMMMMMM               MMMMMMMMEEEEEEEEEEEEEEEEEEEEEE      TTTTTTTTTTT      HHHHHHHHH     HHHHHHHHH     OOOOOOOOO     DDDDDDDDDDDDD <br><br>
</pre>
@par
This is METHOD, a relativistic multi-dimensional, multi-fluid ElectroMagnetoHydroDynamic
solver built and maintained by [Alex Wright](http://cmg.soton.ac.uk/people/ajw1e16/)
under the guidance of [Dr Ian Hawke](https://www.southampton.ac.uk/maths/about/staff/ih3.page).
METHOD is being developed as a result of a PhD project to model neutron star mergers
with multi-fluid models. As a result, ideal and resistive single fluid models exist
that are more conventional in astrophysical model, and a relatively new two-fluid
model adapted from Amano 2016.

@section gettingStarted Getting started
@par
The easiest way to install METHOD is to clone the github [repository](https://github.com/AlexJamesWright/METHOD).
The repository contains all the source files, tests and makefiles that one needs
to get going. To start, copy the HTTPS tag on the main page or run

    git clone https://github.com/AlexJamesWright/METHOD

to copy all files to your machine. To build the correct paths to the data directories
required for saving simulation output, from the METHOD root directory run the
shell script

    bash makePaths.sh

It is then a good idea to run all the tests to ensure the latest build is working. First,
enter the Tests/Serial directory and run

    make test

followed by the same command in Tests/Parallel. The parallel tests will only
run if you have the CUDA drivers installed and an NVIDIA graphics cards. Finally,
from the Tests/Parallel directory, run

    py.test -v Src/comapreParallelAndSerial.py

If all is well, the examples provided will give a good indication of how to use
METHOD.

NOTE: More information about setting up and using METHOD can be found in the project's
`README.md`.
