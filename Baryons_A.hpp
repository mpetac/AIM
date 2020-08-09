#pragma once

#include <complex>
#include "Baryons.hpp"
#include "Parametric_funcs.hpp"
#include "Structs.hpp"

class Baryons_A: public Baryons {

private:
    disk_3p disk1;
    disk_3p disk2;
    bulge_2p bulge;

public:
    Baryons_A(const disk_3p &MN1, const disk_3p &MN2, const bulge_2p &H) {
        disk1 = MN1;
        disk2 = MN2;
        bulge = H;
    }
    
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_d2R2(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2) override;
};
