#include "Halo_gNFW.hpp"

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

Halo_gNFW::Halo_gNFW(const halo_6p& gNFW, const halo_rot_2p& gNFW_rot) {
    if(gNFW.alpha != 1.) std::cout << "Warning: gNFW halo requires alpha = 1!" << std::endl;
    if(gNFW.beta != 3.) std::cout << "Warning: gNFW halo requires beta = 3!" << std::endl;
    if(gNFW.q != 1.) std::cout << "Warning: gNFW halo requires q = 1!" << std::endl;
    Halo_gNFW::halo = gNFW;
    Halo_gNFW::halo_rot = gNFW_rot;
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    double arg = std::arg(r);
    if (arg == 0) r = std::polar(std::abs(r), -1e-8);
    return Numeric_funcs::psi_gNFW(&(Halo_gNFW::halo), r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    double arg = std::arg(r);
    if (arg == 0) r = std::polar(std::abs(r), -1e-8);
//     if (arg > 0) r = std::polar(std::abs(r), arg + M_PI);
    return Numeric_funcs::psi_gNFW_dr2(&(Halo_gNFW::halo), r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    double arg = std::arg(r);
    if (arg == 0) r = std::polar(std::abs(r), -1e-8);
//     if (arg > 0) r = std::polar(std::abs(r), arg + M_PI);
    return Numeric_funcs::psi_gNFW_dr2(&(Halo_gNFW::halo), r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    double arg = std::arg(r);
    if (arg == 0) r = std::polar(std::abs(r), -1e-8);
//     if (arg > 0) r = std::polar(std::abs(r), arg + M_PI);
    return Numeric_funcs::psi_gNFW_d2r2(this, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    double arg = std::arg(r);
    if (arg == 0) r = std::polar(std::abs(r), -1e-8);
//     if (arg > 0) r = std::polar(std::abs(r), arg + M_PI);
    return Numeric_funcs::psi_gNFW_d2r2(this, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / std::pow(Halo_gNFW::halo.q, 2);
    return Parametric_funcs::rho_sABC(m2, Halo_gNFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / std::pow(Halo_gNFW::halo.q, 2);
    return Parametric_funcs::rho_sABC_dz2(m2, Halo_gNFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_gNFW::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / std::pow(Halo_gNFW::halo.q, 2);
    return Parametric_funcs::rho_sABC_d2z2(m2, Halo_gNFW::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 */

std::complex<double> Halo_gNFW::v_phi(std::complex<double> R2) {
    return Parametric_funcs::v_phi(R2, Halo_gNFW::halo_rot);
}


bool Halo_gNFW::is_rotating() {
    if (Halo_gNFW::halo_rot.omega == 0) return 0;
    else return 1;
}
