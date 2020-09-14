#include "Halo_sABC.hpp"

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / Halo_sABC::halo.q2;
    return Parametric_funcs::rho_sABC(m2, Halo_sABC::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / Halo_sABC::halo.q2;
    return Parametric_funcs::rho_sABC_dz2(m2, Halo_sABC::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Halo_sABC::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> m2 = R2 + z2 / Halo_sABC::halo.q2;
    return Parametric_funcs::rho_sABC_d2z2(m2, Halo_sABC::halo);
}

/** 
 * @param R2 Value of the R-coordinate squared
 */

std::complex<double> Halo_sABC::v_phi(std::complex<double> R2) {
    return Parametric_funcs::v_phi(R2, Halo_sABC::halo_rot);
}


bool Halo_sABC::is_rotating() {
    if (Halo_sABC::halo_rot.omega == 0) return 0;
    else return 1;
}
