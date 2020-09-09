#pragma once

#include <complex>
#include "Baryons.hpp"
#include "../Parametric_funcs.hpp"
#include "../Structs.hpp"

/**
 * Implements a baryonic model consisting of two Myiamoto-Nagai disks and a Hernquist bulge.
 */

class Baryons_H_2MN: public Baryons {

private:
    /// Struct with parameters describing the first baryonic disk
    disk_3p disk1;
    /// Struct with parameters describing the second baryonic disk
    disk_3p disk2;
    /// Sturct with parameters describing the bulge
    bulge_2p bulge;

public:
    /// Initializes the baryonic model
    /**
     * @param MN1 Struct with parameters describing the first baryonic disk
     * @param MN2 Struct with parameters describing the second baryonic disk
     * @param H Sturct with parameters describing the bulge
     */
    Baryons_H_2MN(const disk_3p &MN1, const disk_3p &MN2, const bulge_2p &H) {
        disk1 = MN1;
        disk2 = MN2;
        bulge = H;
    }
    
    /// Computes the baryonic gravitational potential
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the first derivative of baryonic gravitational potential with respect to R2
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the first derivative of baryonic gravitational potential with respect to z2
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the second derivative of baryonic gravitational potential with respect to R2 and z2
    std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the second derivative of baryonic gravitational potential with respect to z2
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
};
