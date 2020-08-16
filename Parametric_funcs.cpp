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

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::psi_NFW(std::complex<double> r, const halo_2p& obj) {
    if (r == 0.) return 4. * Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s, 2);
    else return 4. * Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s, 3) * std::log(1. + r / obj.r_s) / r;
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::psi_NFW_dr2(std::complex<double> r, const halo_2p& obj) {
    return 2. * Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s / r, 3) * (r / (r + obj.r_s) - std::log(1. + r / obj.r_s));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::psi_NFW_d2r2(std::complex<double> r, const halo_2p& obj) {
    return Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s, 3) * (3. * std::pow(r + obj.r_s, 2) * std::log(1. + r / obj.r_s) - 4. * std::pow(r, 2) - 3. * r * obj.r_s) / (std::pow(r, 5) * std::pow(r + obj.r_s, 2));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::psi_BUR(std::complex<double> r, const halo_2p& obj) {
    if (r == 0.) return std::pow(Parametric_funcs::PI, 2) * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s, 2);
    else {
        std::complex<double> x = r / obj.r_s;
        return Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s, 2) / x * (2. * (1. + x) * (std::atan(1. / x) + std::log(1. + x)) + (1. - x) * std::log(1. + std::pow(x, 2)) - 1. * Parametric_funcs::PI);
        return 0;
    }
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::psi_BUR_dr2(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return 0.5 * Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * (1. * Parametric_funcs::PI - 2. * std::atan(1. / x) - std::log(1. + std::pow(x, 2)) - 2. * std::log(1. + x)) * std::pow(x, -3);
    return 0;
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::psi_BUR_d2r2(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return 0.25 * Parametric_funcs::PI * Parametric_funcs::G * obj.rho_s * std::pow(obj.r_s, -2) * std::pow(x, -5) / ((1. + x) * (1. + std::pow(x, 2))) * (3. * (1. + x) * (1. + std::pow(x, 2)) * (2. * std::atan(1. / x) + std::log(1. + std::pow(x, 2)) + 2. * std::log(1. + x)) - 3. * Parametric_funcs::PI * (std::pow(x, 3) + std::pow(x, 2) + x + 1.) - 4. * std::pow(x, 3));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::rho_NFW(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return obj.rho_s / (x * std::pow(1. + x, 2));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::rho_NFW_dr2(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return - Parametric_funcs::rho_NFW(r, obj) * std::pow(obj.r_s, -2) * (1. + 3. * x) / (2. * std::pow(x, 2) * (1. + x));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::rho_NFW_d2r2(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return 0.25 * Parametric_funcs::rho_NFW(r, obj) * std::pow(r, -4) * std::pow(1. + x, -2) * 3. * (1. + 4. * x + 5. * std::pow(x, 2));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::rho_BUR(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return obj.rho_s / ((1. + x) * (1. + std::pow(x, 2)));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::rho_BUR_dr2(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return -0.5 * obj.rho_s * std::pow(obj.r_s, -2) * (1. + 2. * x + 3. * std::pow(x, 2)) / (x * std::pow(1. + x + std::pow(x, 2) + std::pow(x, 3), 2));
}

/** 
 * @param r Value of the radial distance
 * @param obj Struct with the halo parameters
 */

std::complex<double> Parametric_funcs::rho_BUR_d2r2(std::complex<double> r, const halo_2p& obj) {
    std::complex<double> x = r / obj.r_s;
    return 0.25 * obj.rho_s * std::pow(obj.r_s, -4) * (1. + 3. * x * (1. + x * (2. + x * (6. + x * (7. + 5. * x))))) * std::pow(x * (1. + x + std::pow(x, 2) + std::pow(x, 3)), -3);
}
