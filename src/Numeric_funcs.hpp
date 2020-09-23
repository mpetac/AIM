#pragma once

#include <iostream>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "Structs.hpp"
#include "halos/Halo.hpp"

/**
 * Class containing the parametric functions (and their relevant derivatives) of the implemented gravitational potentials and density profiles.
 */

class Numeric_funcs {
    
private:
    /// Gravitational constant
    constexpr static double G = 4.3e-6;
    /// Number of intervals used in the numerical integration
    constexpr static size_t nIntervals = 1e5;
    
    /// Utility function for computing the continued fractions used in the evaluation of beta function
    static std::complex<double> beta_cf(std::complex<double> x, double a, double b, int N, int iterations);
    
    /// GSL error handler
    static void GSL_error_func(const char * reason, const char * file, int line, int gsl_errno);

public:
    /// Computes the spheroid's gravitational potential at the origin
    static double psi_spheroid_0(Halo *halo, double q, double tolerance=1e-5);
    
    /// Computes the spheroid's gravitational potential difference with respect to the origin
    static std::complex<double> psi_spheroid(Halo *halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance=1e-3);
    
    /// Computes the first derivative of spheriod's gravitational potential with respect to R2
    static std::complex<double> psi_spheroid_dR2(Halo *halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance=1e-3);
    
    /// Computes the first derivative of spheriod's gravitational potential with respect to z2
    static std::complex<double> psi_spheroid_dz2(Halo *halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance=1e-3);
    
    /// Computes the second derivative of spheriod's gravitational potential with respect to z2
    static std::complex<double> psi_spheroid_d2z2(Halo *halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance=1e-3);
    
    /// Computes the first derivative of spheriod's gravitational potential with respect to R2 and z2
    static std::complex<double> psi_spheroid_d2R2z2(Halo *halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance=1e-3);
    
    /// Computes the gNFW gravitational potential
    static std::complex<double> psi_gNFW(halo_6p *params, std::complex<double> r, int iterations=100);
    
    /// Computes the first derivative of gNFW gravitational potential with respect to r2
    static std::complex<double> psi_gNFW_dr2(halo_6p *params, std::complex<double> r, int iterations=100);
    
    /// Computes the second derivative of gNFW gravitational potential with respect to r2
    static std::complex<double> psi_gNFW_d2r2(Halo *halo, std::complex<double> r);

};
