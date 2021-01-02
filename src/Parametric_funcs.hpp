#pragma once

#include <complex>
#include "Structs.hpp"

/**
 * Class containing the parametric functions (and their relevant derivatives) of the implemented gravitational potentials and density profiles.
 */

class Parametric_funcs {
    
private:
    /// Gravitational constant
    constexpr static double G = 4.3e-6;

public:
    
    /// Computes the Myiamoto-Nagain gravitational potential
    static std::complex<double> psi_MN(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the first derivative of Myiamoto-Nagai gravitational potential with respect to R2
    static std::complex<double> psi_MN_dR2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the first derivative of Myiamoto-Nagai gravitational potential with respect to z2
    static std::complex<double> psi_MN_dz2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the second derivative of Myiamoto-Nagai gravitational potential with respect to z2
    static std::complex<double> psi_MN_d2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the first derivative of Myiamoto-Nagai gravitational potential with respect to R2 and z2
    static std::complex<double> psi_MN_d2R2z2(std::complex<double> R2, std::complex<double> z2, const struct disk_3p& obj);
    
    /// Computes the Hernquist gravitational potential
    static std::complex<double> psi_H(std::complex<double> r, const struct bulge_2p& obj);
    
    /// Computes the first derivative of Hernquist gravitational potential with respect to r2
    static std::complex<double> psi_H_dr2(std::complex<double> r, const struct bulge_2p& obj);
    
    /// Computes the second derivative of Hernquist gravitational potential with respect to r2
    static std::complex<double> psi_H_d2r2(std::complex<double> r, const struct bulge_2p& obj);
    
    /// Computes the NFW gravitational potential
    static std::complex<double> psi_NFW(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the first derivative of NFW gravitational potential with respect to r2
    static std::complex<double> psi_NFW_dr2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the second derivative of NFW gravitational potential with respect to r2
    static std::complex<double> psi_NFW_d2r2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the Burkert gravitational potential
    static std::complex<double> psi_BUR(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the first derivative of Burkert gravitational potential with respect to r2
    static std::complex<double> psi_BUR_dr2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the second derivative of Burkert gravitational potential with respect to r2
    static std::complex<double> psi_BUR_d2r2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the NFW density profile
    static std::complex<double> rho_NFW(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the first derivative of NFW density profile with respect to r2
    static std::complex<double> rho_NFW_dr2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the second derivative of NFW density profile with respect to r2
    static std::complex<double> rho_NFW_d2r2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the Burkert density profile
    static std::complex<double> rho_BUR(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the first derivative of Burkert density profile with respect to r2
    static std::complex<double> rho_BUR_dr2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the second derivative of Burkert density profile with respect to r2
    static std::complex<double> rho_BUR_d2r2(std::complex<double> r, const struct halo_2p& obj);
    
    /// Computes the alpha-beta-gamma density profile
    static std::complex<double> rho_sABC(std::complex<double> m2, const struct halo_6p& obj);
    
    /// Computes the first derivative of alpha-beta-gamma density profile with respect to r2
    static std::complex<double> rho_sABC_dz2(std::complex<double> m2, const struct halo_6p& obj);
    
    /// Computes the second derivative of alpha-beta-gamma density profile with respect to r2
    static std::complex<double> rho_sABC_d2z2(std::complex<double> m2, const struct halo_6p& obj);
    
    /// Computes the rotational velocity
    static std::complex<double> v_phi(std::complex<double> R2, std::complex<double> z2, const struct halo_rot_3p& obj);
    
    /// Computes the first derivative of rotational velocity with respect to z2
    static std::complex<double> v_phi_dz2(std::complex<double> R2, std::complex<double> z2, const struct halo_rot_3p& obj);
    
    /// Computes the second derivative of rotational velocity with respect to z2
    static std::complex<double> v_phi_d2z2(std::complex<double> R2, std::complex<double> z2, const struct halo_rot_3p& obj);
};
