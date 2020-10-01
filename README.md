# AIM - axisymmetric inversion method

Code for computing equilibrium axisymmetric phase-space distribution function. It provides highly efficent and user-friendly implementation of an algorithm presented in [Hunter \& Qian](https://academic.oup.com/mnras/article/262/2/401/1161204), [Qian et al.](https://academic.oup.com/mnras/article/274/2/602/2896126) and [Petac \& Ullio](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.043003). The code is written in C++, however, it can be also be deployed from Python using the provided pybind11 wrapper.

## Requirements

To compile the source code the following tools have to be installed on your system:
* g++ (recommended version 7.5)
* GSL (recommended version 2.4)
* cmake 2.6 or later
* pybind11 2.5 or later (optional)

## Usage

The code can be used either from C++ or Python. While the Python wrapper provides modern and simple interface, running it from C++ allows for more options (mainly related to the accuracy of various numerical perscriptions).
