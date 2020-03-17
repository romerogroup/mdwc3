Molecular Dynamics With Constraints (mdwc3) - Python 3 version
=============================================================

**This is the Python 3 version of the original mdwc code**

The molecular dynamics with constraints (mdwc3) package is a command line open source Python program for molecular dynamics simulations. It does constraint molecular dynamics with either NPT (keeping pressure constant with the Parrinello Rahman lagrangian, and keeping the temperature constant with the Nose thermostat) or NVT (keeping the temperature constant with the Nose thermostat). mdwc3 performs constraint dynamics following the SHAKE algorithm and allows constraints for bond distances, angles, atomic positions, lattice parameters (a, b, c), angles between lattice vectors, and volume of the unit cell.

Mdwc has the mdwc3 script to interface the mdwc code with an external first principles code for the calculations of energies and forces, so far, the only supported DFT code is [Abinit](https://www.abinit.org/). You can use the mdwc functions to carry on constrained molecular dynamics with other first principals codes (for more details look at the jupyter notebook section)


mdwc requirements
-----------------
For its correct functioning mlmd needs the following codes:

1. [Python](https://www.python.org/download/releases/3.7/ "Python") = 3.7 

2. [Numpy](http://www.numpy.org/ "Numpy") >= 1.17.2

* [Abinit](https://www.abinit.org/) to use the mdwc3 script

mdwc user manual
----------------

This is for the Python 2 version, but it is applicable to the Python 3 version as well. 

[mdwc](https://molecular-dynamics-with-constraints.github.io/)

installation
------------

First, download or clone the repository.

    git clone https://github.com/romerogroup/mdwc3.git

Then run the setup.

    python setup.py install

Or simply install it with pip.

    pip install mdwc3

Running mdwc3
------------

    mdwc3 -name c_8 -mpirun True -np 32 -NVT False    

Contributors
------------
* Prof. Aldo H. Romero [West Virginia University] 

* Arturo Hernandez [West Virginia University] (Developer)

* Uthpala Herath   [West Virginia University] (Simulation and testing) 

* Pedram Tavazohi  [West Virginia University] (Simulation and testing)
