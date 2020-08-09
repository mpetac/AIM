#include "Parametric_funcs.hpp"

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param obj Struct with the Miyamoto-Nagai parameters
 */

std::complex<double> Parametric_funcs::psi_MN(std::complex<double> R2, std::complex<double> z2, const disk_3p& obj) {
    return Parametric_funcs::G * obj.M / std::sqrt(R2 + std::pow(obj.a + std::sqrt(std::pow(obj.b, 2) + z2), 2));
}

/**
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param obj Struct with the Miyamoto-Nagai parameters
 */

std::complex<double> Parametric_funcs::psi_MN_dR2(std::complex<double> R2, std::complex<double> z2, const disk_3p& obj) {
    return - Parametric_funcs::G * obj.M / (2. * std::pow(R2 + std::pow(obj.a + std::sqrt(std::pow(obj.b, 2) + z2), 2), 1.5));
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param obj Struct with the Miyamoto-Nagai parameters
 */

std::complex<double> Parametric_funcs::psi_MN_dz2(std::complex<double> R2, std::complex<double> z2, const disk_3p& obj) {
    std::complex<double> sqrt_b2z2 = std::sqrt(std::pow(obj.b, 2) + z2);
    return - Parametric_funcs::G * obj.M * (obj.a + sqrt_b2z2) / (2. * sqrt_b2z2 * std::pow(R2 + std::pow(obj.a + sqrt_b2z2, 2), 1.5));
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param obj Struct with the Miyamoto-Nagai parameters
 */

std::complex<double> Parametric_funcs::psi_MN_d2z2(std::complex<double> R2, std::complex<double> z2, const disk_3p& obj) {
    std::complex<double> sqrt_b2z2 = std::sqrt(std::pow(obj.b, 2) + z2);
    std::complex<double> csqrt_b2z2 = std::pow(std::pow(obj.b, 2) + z2, 1.5);
    return Parametric_funcs::G * obj.M * (std::pow(obj.a, 3) + 5. * std::pow(obj.a , 2) * sqrt_b2z2 + 3. * csqrt_b2z2 + obj.a * (7 * std::pow(obj.b, 2) + R2 + 7. * z2)) / (4. * csqrt_b2z2 * std::pow(R2 + std::pow(obj.a + sqrt_b2z2, 2), 2.5));
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param obj Struct with the Miyamoto-Nagai parameters
 */

std::complex<double> Parametric_funcs::psi_MN_d2R2z2(std::complex<double> R2, std::complex<double> z2, const disk_3p& obj) {
    std::complex<double> sqrt_b2z2 = std::sqrt(std::pow(obj.b, 2) + z2);
    return 3. * Parametric_funcs::G * obj.M * (obj.a + sqrt_b2z2) / (4. * sqrt_b2z2 * std::pow(R2 + std::pow(obj.a + sqrt_b2z2, 2), 2.5));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the Herenquist parameters
 */

std::complex<double> Parametric_funcs::psi_H(std::complex<double> r, const bulge_2p& obj) {
    return Parametric_funcs::G * obj.M / (r + obj.a);
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the Herenquist parameters
 */

std::complex<double> Parametric_funcs::psi_H_dr2(std::complex<double> r, const bulge_2p& obj) {
    return - Parametric_funcs::G * obj.M / (2. * r * std::pow(r + obj.a, 2));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the Herenquist parameters
 */

std::complex<double> Parametric_funcs::psi_H_d2r2(std::complex<double> r, const bulge_2p& obj) {
    return Parametric_funcs::G * obj.M * (3. * r + obj.a) / (4. * std::pow(r, 3) * std::pow(r + obj.a, 3));
}
