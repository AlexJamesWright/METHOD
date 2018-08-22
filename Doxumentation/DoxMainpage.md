
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
  A few example simulations have been provided that illustrate how to use the
various classes. By typing `make run` in one of the example directories, the
relevant object files will be built and executed. Data is saved in the *Examples/Data*
directory and is easily viewed using the interactivePlot script, run from the
root Example directory with something like `spyder interactivePlot.py`. For the
Kelvin-Helmholtz simulation, running the `animation.py` script will create an
animatation called `Output.gif` in the root Example directory to view (may take up
to ten minutes to run the simulation and make the animation).
  @note When generating animations, besure to delete all TimeSeries data after each run.
