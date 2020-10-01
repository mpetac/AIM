# AIM - axisymmetric inversion method

Code for computing equilibrium axisymmetric phase-space distribution function. It provides highly efficent and user-friendly implementation of an algorithm presented in [Hunter \& Qian](https://academic.oup.com/mnras/article/262/2/401/1161204), [Qian et al.](https://academic.oup.com/mnras/article/274/2/602/2896126) and [Petac \& Ullio](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.043003) (please cite the latter work if you use this code in your work). The code is written in C++, however, it can be also be deployed from Python using the provided pybind11 wrapper.

## Requirements

To compile the source code the following tools have to be installed on your system:
* g++ (recommended version 7.5)
* GSL (recommended version 2.4)
* cmake 2.6 or later
* pybind11 2.5 or later (optional, needed for compiling the Python wrapper)
* doxygen (optional, needed for automatic generation of the documentation)

## Usage

The code can be used either from C++ or Python. While the Python wrapper provides modern and simple interface, running it from C++ allows for more control (mainly related to the accuracy of various numerical perscriptions).

#### C++

The code can be used in C++ as demonstrated in the following example (see main.cpp):


```C++
    //Include the relevant 
    #include "src/Observables.hpp"
    #include "src/Inversion.hpp"
    #include "src/Model.hpp"
    #include "src/Structs.hpp"
    #include "src/halos/Halo_NFW.hpp"
    #include "src/halos/Halo_gNFW.hpp"
    #include "src/halos/Halo_sABC.hpp"
    #include "src/baryons/Baryons_H_2MN.hpp"


    // Define struct with DM halo parameters. Here we asume spherical DM density profile with density 1e7 M_sol / kpc^3 and scale density of 13 kpc.
    halo_2p p_nfw = {1e7, 13.};
    // Initialize the DM halo object.
    Halo_NFW halo(p_nfw);

    // Define structs related to the baryonic distribution. In this example we assume a model consisting of tow Myiamoto-Nagai disks and a spherical Hernquist bulge.
    disk_3p disk1 = {5e10, 3.6, 0.5};
    disk_3p disk2 = {0., 1., 1.};
    bulge_2p bulge = {1e11, 1.};
    // Initialize the baryonic model
    Baryons_H_2MN baryons(disk1, disk2, bulge);

    // Initialize the galactic model using the previously defined halo and baryons
    Model model(&halo, &baryons);

    // Interpolate the PSDF obtained for the specified galactic model with given number of relative energy and angular momentum points
    Inversion psdf(&model, 1000, 10);

    // Initialize the class for computing various observable quantities from the PSDF (namely DM density and various projections of the velocity distribution)
    Observables obs(&model, &psdf);
```
