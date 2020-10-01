# AIM - axisymmetric inversion method

Code for computing equilibrium axisymmetric phase-space distribution function. It provides highly efficent and user-friendly implementation of an algorithm presented in [Hunter \& Qian](https://academic.oup.com/mnras/article/262/2/401/1161204), [Qian et al.](https://academic.oup.com/mnras/article/274/2/602/2896126) and [Petac \& Ullio](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.043003). The code is written in C++, however, it can be also be deployed from Python using the provided wrapper module.

This project is distributed under GPL v3.0 licence. If you use this code for a scientific project, please cite our paper that discusses the details of the implemented inversion method: [Petac \& Ullio](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.043003).

## Requirements

To compile the source code the following tools have to be installed on your system:
* gcc (recommended version 7.5 or later)
* GSL (recommended version 2.4 or later)
* cmake version 2.6 or later
* pybind11 version 2.5 or later (optional, needed for compiling the Python wrapper)
* doxygen (optional, needed for automatic generation of the documentation)

## Examples

The code can be used either from C++ or Python. While the Python wrapper provides modern and simple interface, running it from C++ allows for more control (mainly related to the accuracy of various numerical perscriptions). Futhermore, the C++ class structure allows for simple extension the code to new DM density profile, halo rotational properties and baryonic gravitational potentials.

#### C++

The code can be used in C++ as demonstrated in the following example (see main.cpp):


```c++
//Include the relevant source files.
#include "src/Observables.hpp"
#include "src/Inversion.hpp"
#include "src/Model.hpp"
#include "src/Structs.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/halos/Halo_gNFW.hpp"
#include "src/halos/Halo_sABC.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"


// Define struct with DM halo parameters.
// Here we asume spherical NFW density profile.
// First parameter is the scale density (in units of M_sol / kpc^3) while the second parameter is scale density (in units of kpc).
halo_2p p_nfw = {1e7, 13.};
// Initialize the DM halo object.
Halo_NFW halo(p_nfw);

// Define structs related to the baryonic distribution.
// In this example we assume a model consisting of two Myiamoto-Nagai disks (first parameter is the disk mass in units of M_sol,+ while the second and third parameters are scale length and scale height in units of kpc) and a spherical Hernquist bulge (first parameter is the bulge mass in units of M_sol while the second parameter is the scale lenght in units of kpc).
disk_3p disk1 = {5e10, 3.6, 0.5};
disk_3p disk2 = {0., 1., 1.};
bulge_2p bulge = {1e11, 1.};
// Initialize the baryonic model
Baryons_H_2MN baryons(disk1, disk2, bulge);

// Initialize the galactic model using the previously defined halo and baryons.
Model model(&halo, &baryons);

// Interpolate the PSDF obtained for the specified galactic model with given number of relative energy and angular momentum points.
Inversion psdf(&model, 1000, 10);

// Initialize the class for computing various observable quantities from the PSDF (namely DM density and various projections of the velocity distribution).
Observables obs(&model, &psdf);

// Evaluate the DM density (computed numerically from the obtained distribution function) in 20 points along the radial direciton.
int nPts = 20;
double logRmin = -1, logRmax = 3.;
double Rpts[nPts], zpts[nPts], result[nPts];
for (int i = 0; i < nPts; i++) {
    Rpts[i] = std::pow(10., logRmin + (logRmax - logRmin) * i / (nPts - 1));
    zpts[i] = 0;
}
obs.rho(nPts, Rpts, zpts, result);

// Evaluate the probability distribution for the velocity magnitude at R=8.122 kpc and z=0 using 100 velocity points.
int nVel = 100;
double pv_mag[2 * nVel];
obs.pv_mag(nVel, 8.122, 0, pv_mag);
```

To add new baryonic or halo models add the corresponding source files in /src/halo or /src/baryons folders, respectfully. They have to extend the virtual parent class Halo.hpp or Baryons.hpp.

#### Python

To use the Python wrapper simply copy the AIM.so file in your project's directory. Example of its usage:

```python
import numpy as np
# Import the wrapper module.
import AIM

# Initialize a new instance of the phase-space distribution function (the argument controls weather you want verbose output or not).
f = AIM.PSDF(True)

# Setup a halo model with spherical NFW density profile.
f.setHalo('NFW', [1e7, 13.])
# Setup a baryonic model composed of Hernquist bulge and two Myiamoto-Nagai disks.
f.setBaryons('H_2MN', [5e10, 3.6, 0.5, 1e10, 2.5, 1., 1e10, 0.5])

# Compute the distribution function with a given number of relative energy and angular momentum points.
f.compute(1000, 20)

# Compute the DM density profile from the obtained distribution function in 20 radial points.
Rpts = np.linspace(1, 200, 20)
zpts = np.zeros(20)
density_profile = f.rho(20, Rpts, zpts)

# Compute the local (R~8.122 kpc, z~0) velocity distribution of DM using 100 velocity points.
velocity_distribution = f.pv_mag(100, 8.122, 0)
```

## Documentation

Source files contain brief comments that should be sufficient for understanding of the code. Comprehensive documentation based on these comments can be generated through Doxygen using the provided setup file ("Doxyfile").
