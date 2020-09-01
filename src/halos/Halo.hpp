#pragma once

#include <complex>

class Halo {
    
public:
    virtual std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    virtual std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    virtual std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    virtual std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    virtual std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    virtual std::complex<double> rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    virtual std::complex<double> rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    virtual std::complex<double> rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) = 0;
    
    virtual std::complex<double> v_phi(std::complex<double> R2) = 0;
    
    virtual bool is_rotating() = 0;
};
