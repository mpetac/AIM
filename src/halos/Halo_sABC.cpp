#include "Halo_sABC.hpp"

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

Halo_sABC::Halo_sABC(const halo_6p& sABC, const halo_rot_3p& sABC_rot) {
    Halo_sABC::halo = sABC;
    Halo_sABC::halo_rot = sABC_rot;
    Halo_sABC::psi0 = std::real(Numeric_funcs::psi_spheroid_0(this, Halo_sABC::halo.q));
}


std::complex<double> Halo_sABC::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    if (R2 == 0. && z2 == 0. && r == 0.) return Halo_sABC::psi0;
    return Halo_sABC::psi0 - Numeric_funcs::psi_spheroid(this, R2, z2, Halo_sABC::halo.q);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    if (R2 == 0. && z2 == 0. && r == 0.) return 0.;
    return Numeric_funcs::psi_spheroid_dR2(this, R2, z2, Halo_sABC::halo.q);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    if (R2 == 0. && z2 == 0. && r == 0.) return 0.;
    return Numeric_funcs::psi_spheroid_dz2(this, R2, z2, Halo_sABC::halo.q);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    if (R2 == 0. && z2 == 0. && r == 0.) return 0.;
    return Numeric_funcs::psi_spheroid_d2R2z2(this, R2, z2, Halo_sABC::halo.q);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    if (R2 == 0. && z2 == 0. && r == 0.) return 0.;
    return Numeric_funcs::psi_spheroid_d2z2(this, R2, z2, Halo_sABC::halo.q);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / std::pow(Halo_sABC::halo.q, 2);
    return Parametric_funcs::rho_sABC(m2, Halo_sABC::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / std::pow(Halo_sABC::halo.q, 2);
    return Parametric_funcs::rho_sABC_dz2(m2, Halo_sABC::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / std::pow(Halo_sABC::halo.q, 2);
    return Parametric_funcs::rho_sABC_d2z2(m2, Halo_sABC::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 */

std::complex<double> Halo_sABC::v_phi(std::complex<double> R2, std::complex<double> z2) {
    return Parametric_funcs::v_phi(R2, z2, Halo_sABC::halo_rot);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 */

std::complex<double> Halo_sABC::v_phi_dz2(std::complex<double> R2, std::complex<double> z2) {
    return Parametric_funcs::v_phi_dz2(R2, z2, Halo_sABC::halo_rot);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 */

std::complex<double> Halo_sABC::v_phi_d2z2(std::complex<double> R2, std::complex<double> z2) {
    return Parametric_funcs::v_phi_d2z2(R2, z2, Halo_sABC::halo_rot);
}


bool Halo_sABC::is_rotating() {
    if (Halo_sABC::halo_rot.omega == 0) return 0;
    else return 1;
}
