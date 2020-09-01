#pragma once

#include <complex>
#include "Halo.hpp"
#include "../Parametric_funcs.hpp"
#include "../Structs.hpp"

class Halo_NFW: public Halo {

private:
    halo_2p halo;
    halo_rot_2p halo_rot;

public:
    Halo_NFW(const halo_2p &NFW, const halo_rot_2p &NFW_rot = {0, 0}) {
        halo = NFW;
        halo_rot = NFW_rot;
    }
    
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    std::complex<double> rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    std::complex<double> rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    std::complex<double> rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    std::complex<double> v_phi(std::complex<double> R2) override;
    
    bool is_rotating() override;
};
