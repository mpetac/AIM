#pragma once

#include <complex>

/**
 * Abstract class for describing the galactic baryons.
 */

class Baryons {
    
public:
    /// Computes the baryonic gravitational potential
    virtual std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the first derivative of baryonic gravitational potential with respect to R2
    virtual std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the first derivative of baryonic gravitational potential with respect to z2
    virtual std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the second derivative of baryonic gravitational potential with respect to R2 and z2
    virtual std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the second derivative of baryonic gravitational potential with respect to z2
    virtual std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
};
