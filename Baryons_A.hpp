#pragma once

#include <complex>
#include "Baryons.hpp"
#include "Structs.hpp"

class Baryons_A: public Baryons {

private:
    disk_3p MN1;
    disk_3p MN2;
    bulge_2p H;

public:
    Baryons_A(const disk_3p &disk1, const disk_3p &disk2, const bulge_2p &bulge) {
        MN1 = disk1;
        MN2 = disk2;
        H = bulge;
    }
    
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_d2R2(std::complex<double> R2, std::complex<double> z2) override;
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2) override;
    
    std::complex<double> psi_MN(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    std::complex<double> psi_MN_dR2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    std::complex<double> psi_MN_dz2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    std::complex<double> psi_MN_d2R2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    std::complex<double> psi_MN_d2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    
    /*
    std::complex<double> psi_H(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    std::complex<double> psi_H_dr2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    std::complex<double> psi_H_d2r2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p&);
    */
};
