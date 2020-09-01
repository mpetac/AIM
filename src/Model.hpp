#pragma once

#include <iostream>
#include <complex>
#include <gsl/gsl_multimin.h>
#include "Structs.hpp"
#include "halos/Halo.hpp"
#include "baryons/Baryons.hpp"

class Model {

private:
    bool verbose;

public:
    Halo* halo;
    Baryons* baryons;
    
    Model(Halo* halo_model, Baryons* baryons_model, bool v = 1) {
        halo = halo_model;
        baryons = baryons_model;
        verbose = v;
    }
    
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    std::complex<double> rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    std::complex<double> rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    std::complex<double> rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    std::complex<double> rho_d2psi2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r);
    
    std::complex<double> v_phi(std::complex<double> R2);
    
    std::complex<double> psi_inverse(std::complex<double> xi, double E, double Lz, std::complex<double> z0 = (1., 0), double tolerance = 1e-6, int limit = 200);
    
    double Rcirc (double E, double tolerance = 1e-6, int limit = 1000);
    
    bool is_rotating();
};
