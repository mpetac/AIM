#pragma once

#include <complex>

class Baryons {
    
public:
    virtual std::complex<double> psi(std::complex<double> R2, std::complex<double> z2) = 0;
    virtual std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2) = 0;
    virtual std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2) = 0;
    virtual std::complex<double> psi_d2R2(std::complex<double> R2, std::complex<double> z2) = 0;
    virtual std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2) = 0;
};
