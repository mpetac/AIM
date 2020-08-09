#include "Baryons_A.hpp"

std::complex<double> Baryons_A::psi(std::complex<double> R2, std::complex<double> z2) {
    std::complex<double> r = std::sqrt(R2 + z2);
    return Parametric_funcs::psi_MN(R2, z2, Baryons_A::disk1) + Parametric_funcs::psi_MN(R2, z2, Baryons_A::disk2) + Parametric_funcs::psi_H(r, Baryons_A::bulge);
}

std::complex<double> Baryons_A::psi_dR2(std::complex<double> R2, std::complex<double> z2) {
}

std::complex<double> Baryons_A::psi_dz2(std::complex<double> R2, std::complex<double> z2) {
}

std::complex<double> Baryons_A::psi_d2R2(std::complex<double> R2, std::complex<double> z2) {
}

std::complex<double> Baryons_A::psi_d2z2(std::complex<double> R2, std::complex<double> z2) {
}
