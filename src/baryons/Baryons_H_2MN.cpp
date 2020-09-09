#include "Baryons_H_2MN.hpp"

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Baryons_H_2MN::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_MN(R2, z2, Baryons_H_2MN::disk1) + Parametric_funcs::psi_MN(R2, z2, Baryons_H_2MN::disk2) + Parametric_funcs::psi_H(r, Baryons_H_2MN::bulge);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Baryons_H_2MN::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_MN_dR2(R2, z2, Baryons_H_2MN::disk1) + Parametric_funcs::psi_MN_dR2(R2, z2, Baryons_H_2MN::disk2) + Parametric_funcs::psi_H_dr2(r, Baryons_H_2MN::bulge);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Baryons_H_2MN::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_MN_dz2(R2, z2, Baryons_H_2MN::disk1) + Parametric_funcs::psi_MN_dz2(R2, z2, Baryons_H_2MN::disk2) + Parametric_funcs::psi_H_dr2(r, Baryons_H_2MN::bulge);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Baryons_H_2MN::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_MN_d2R2z2(R2, z2, Baryons_H_2MN::disk1) + Parametric_funcs::psi_MN_d2R2z2(R2, z2, Baryons_H_2MN::disk2) + Parametric_funcs::psi_H_d2r2(r, Baryons_H_2MN::bulge);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of the radial distance
 */

std::complex<double> Baryons_H_2MN::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_MN_d2z2(R2, z2, Baryons_H_2MN::disk1) + Parametric_funcs::psi_MN_d2z2(R2, z2, Baryons_H_2MN::disk2) + Parametric_funcs::psi_H_d2r2(r, Baryons_H_2MN::bulge);
}
