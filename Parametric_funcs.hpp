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
    
    /// Computes the Myiamoto-Nagain potential
    static std::complex<double> psi_MN(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    
    /// Computes the first derivative of Myamoto-Nagai potential with respect to R2
    static std::complex<double> psi_MN_dR2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    
    /// Computes the first derivative of Myamoto-Nagai potential with respect to z2
    static std::complex<double> psi_MN_dz2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    
    /// Computes the second derivative of Myamoto-Nagai potential with respect to z2
    static std::complex<double> psi_MN_d2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    
    /// Computes the first derivative of Myamoto-Nagai potential with respect to R2 and z2
    static std::complex<double> psi_MN_d2R2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    
    /// Computes the Hernquist potential
    static std::complex<double> psi_H(std::complex<double> r, const struct bulge_2p&);
    
    /// Computes the first derivative of Hernquist potential with respect to r2
    static std::complex<double> psi_H_dr2(std::complex<double> r, const struct bulge_2p&);
    
    /// Computes the second derivative of Hernquist potential with respect to r2
    static std::complex<double> psi_H_d2r2(std::complex<double> r, const struct bulge_2p&);
    
};
