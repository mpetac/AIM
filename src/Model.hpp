#pragma once

#include <iostream>
#include <complex>
#include <gsl/gsl_multimin.h>
#include "Structs.hpp"
#include "halos/Halo.hpp"
#include "baryons/Baryons.hpp"

/**
 * Top-level class used for describing the system (i.e. galaxy).
 */

class Model {

private:
    /// Parameter controling verbose output
    bool verbose;

public:
    /// Class describing the DM halo
    Halo* halo;
    /// Class desribing the baryons
    Baryons* baryons;
    /// Normalization of the gravitational potential
    double psi0;
    
    /// Initializer for th galactic model
    Model(Halo* halo_model, Baryons* baryons_model, bool v = 1);
    
    /// Computes the total gravitational potential
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the first derivative of total gravitational potential with respect to R2
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the first derivative of total gravitational potential with respect to z2
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the second derivative of total gravitational potential with respect to R2 and z2
    std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the second derivative of total gravitational potential with respect to z2
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the DM density profile
    std::complex<double> rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the first derivative of DM density profile with respect to z2
    std::complex<double> rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the second derivative of DM density profile with respect to z2
    std::complex<double> rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the second derivative of DM density profile with respect to the gravitational potential
    std::complex<double> rho_d2psi2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    /// Computes the DM rotational velocity
    std::complex<double> v_phi(std::complex<double> R2);
    
    /// Compute the inverse of the total gravitational potential (i.e. z2 = psi^{-1}(xi, R2))
    std::complex<double> psi_inverse(std::complex<double> xi, double E, double Lz, std::complex<double> z0 = 1e3+1e3i, double tolerance = 1e-8, int limit = 200);
    
    /// Computes the circular radius corresponding to a given energy
    double Rcirc (double E, double tolerance = 1e-8, int limit = 1000);
    
    /// Check whether the halo is rotating or not
    bool is_rotating();
};
