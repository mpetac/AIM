#include "Halo_NFW.hpp"

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_dr2(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_dr2(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_d2r2(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_d2r2(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::rho_NFW(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::rho_NFW_dr2(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_NFW::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::rho_NFW_d2r2(r, Halo_NFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 */

std::complex<double> Halo_NFW::v_phi(std::complex<double> R2) {
    return Parametric_funcs::v_phi(R2, Halo_NFW::halo_rot);
}


bool Halo_NFW::is_rotating() {
    if (Halo_NFW::halo_rot.omega == 0) return 0;
    else return 1;
}
