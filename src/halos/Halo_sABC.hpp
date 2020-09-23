#pragma once

#include <complex>
#include "Halo.hpp"
#include "../Numeric_funcs.hpp"
#include "../Parametric_funcs.hpp"
#include "../Structs.hpp"

/**
 * Class implementing the spheroidal alpha-beta-gamma halo model. The halo can be made to rotate "on cylinders" if the appropriate rotational parameters are provided.
 */

class Halo_sABC: public Halo {

private:
    /// Struct with halo parameters
    halo_6p halo;
    /// Struct with halo parameters
    halo_rot_2p halo_rot;
    /// Value of the potential at the center of halo
    double psi0;
    /// Function for computing the value of the potential at the center of the halo
    double psi_0();

public:
    /// Initializes halo with spheroidal alpha-beta-gamma density profile
    /**
     * @param BUR Struct with DM density parameters
     * @param BUR_rot Struct with parameters describing the halo rotation
     */
    Halo_sABC(const halo_6p &sABC, const halo_rot_2p &sABC_rot = {0, 0});
    
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
    std::complex<double> v_phi(std::complex<double> R2) override;
    
    /// Checks if halo is rotating (i.e. has non-zero mean rotation around the z-axis)
    bool is_rotating() override;
};
