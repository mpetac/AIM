#pragma once

#include <complex>

/**
 * Abstract class for describing the DM halo.
 */

class Halo {
    
public:
    /// Computes the halo's gravitational potential
    virtual std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the first derivative of halo's gravitational potential with respect to R2
    virtual std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the first derivative of halo's gravitational potential with respect to z2
    virtual std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the second derivative of halo's gravitational potential with respect to R2 and z2
    virtual std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the second derivative of halo's gravitational potential with respect to z2
    virtual std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the DM density
    virtual std::complex<double> rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the first derivative of DM density with respect to z2
    virtual std::complex<double> rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the second derivative of DM density with respect to z2
    virtual std::complex<double> rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    /// Computes the halo's mean rotation around the z-axis
    virtual std::complex<double> v_phi(std::complex<double> R2) = 0;
    
    /// Checks if halo is rotating (i.e. has non-zero mean rotation around the z-axis)
    virtual bool is_rotating() = 0;
};
