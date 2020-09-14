#pragma once

#include <complex>
#include "Structs.hpp"

/**
 * Class containing the parametric functions (and their relevant derivatives) of the implemented gravitational potentials and density profiles.
 */

class Parametric_funcs {
    
private:
    /// Gravitational constant
    constexpr static double G = 4.3e-6;

public:
    
    /// Computes the Myiamoto-Nagain gravitational potential
    static std::complex<double> psi_spheroid(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the first derivative of Myiamoto-Nagai gravitational potential with respect to R2
    static std::complex<double> psi_spheroid_dR2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the first derivative of Myiamoto-Nagai gravitational potential with respect to z2
    static std::complex<double> psi_spheroid_dz2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the second derivative of Myiamoto-Nagai gravitational potential with respect to z2
    static std::complex<double> psi_spheroid_d2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the first derivative of Myiamoto-Nagai gravitational potential with respect to R2 and z2
    static std::complex<double> psi_spheroid_d2R2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
};
