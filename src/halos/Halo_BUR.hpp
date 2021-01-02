#pragma once

#include <complex>
#include "Halo.hpp"
#include "../Parametric_funcs.hpp"
#include "../Structs.hpp"

/**
 * Class implementing the spherical Burkert halo. The halo can be made to rotate "on cylinders" if the appropriate rotational parameters are provided.
 */

class Halo_BUR: public Halo {

private:
    /// Struct with halo parameters
    halo_2p halo;
    /// Struct with halo rotation parameters
    halo_rot_3p halo_rot;

public:
    /// Initializes halo with Burkert density profile
    /**
     * @param BUR Struct with DM density parameters
     * @param BUR_rot Struct with parameters describing the halo rotation
     */
    Halo_BUR(const halo_2p &BUR, const halo_rot_3p &BUR_rot = {0, 1., 1.}) {
        halo = BUR;
        halo_rot = BUR_rot;
    }
    
    /// Computes the halo's gravitational potential
    std::complex<double> psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the first derivative of halo's gravitational potential with respect to R2
    std::complex<double> psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the first derivative of halo's gravitational potential with respect to z2
    std::complex<double> psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the second derivative of halo's gravitational potential with respect to R2 and z2
    std::complex<double> psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the second derivative of halo's gravitational potential with respect to z2
    std::complex<double> psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the DM density
    std::complex<double> rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the first derivative of DM density with respect to z2
    std::complex<double> rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the second derivative of DM density with respect to z2
    std::complex<double> rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) override;
    
    /// Computes the halo's mean rotation around the z-axis
    std::complex<double> v_phi(std::complex<double> R2, std::complex<double> z2) override;
    
    /// Computes the first derivative of halo's mean rotation around the z-axis with respect to z2
    std::complex<double> v_phi_dz2(std::complex<double> R2, std::complex<double> z2) override;
    
    /// Computes the second derivative of halo's mean rotation around the z-axis with respect to z2
    std::complex<double> v_phi_d2z2(std::complex<double> R2, std::complex<double> z2) override;
    
    /// Checks if halo is rotating (i.e. has non-zero mean rotation around the z-axis)
    bool is_rotating() override;
};
