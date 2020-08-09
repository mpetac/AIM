#include <complex>
#include "Baryons_A.hpp"
#include "Structs.hpp"

std::complex<double> Baryons_A::psi(std::complex<double> R2, std::complex<double> z2) {
    return psi_MN(R2, z2, Baryons_A::MN1) + psi_MN(R2, z2, Baryons_A::MN2);
}

std::complex<double> Baryons_A::psi_dR2(std::complex<double> R2, std::complex<double> z2) {
    return psi_MN(R2, z2, Baryons_A::MN1) + psi_MN(R2, z2, Baryons_A::MN2);
}

std::complex<double> Baryons_A::psi_dz2(std::complex<double> R2, std::complex<double> z2) {
    return psi_MN(R2, z2, Baryons_A::MN1) + psi_MN(R2, z2, Baryons_A::MN2);
}

std::complex<double> Baryons_A::psi_d2R2(std::complex<double> R2, std::complex<double> z2) {
    return psi_MN(R2, z2, Baryons_A::MN1) + psi_MN(R2, z2, Baryons_A::MN2);
}

std::complex<double> Baryons_A::psi_d2z2(std::complex<double> R2, std::complex<double> z2) {
    return psi_MN(R2, z2, Baryons_A::MN1) + psi_MN(R2, z2, Baryons_A::MN2);
}

std::complex<double> Baryons_A::psi_MN(std::complex<double> R2, std::complex<double> z2, const disk_3p& disk) {
    return disk.M / std::sqrt(R2 + std::pow(disk.a + std::sqrt(std::pow(disk.b, 2) + z2), 2));
}

std::complex<double> Baryons_A::psi_MN_dR2(std::complex<double> R2, std::complex<double> z2, const disk_3p& disk) {
    return disk.M / std::sqrt(R2 + std::pow(disk.a + std::sqrt(std::pow(disk.b, 2) + z2), 2));
}

std::complex<double> Baryons_A::psi_MN_dz2(std::complex<double> R2, std::complex<double> z2, const disk_3p& disk) {
    return disk.M / std::sqrt(R2 + std::pow(disk.a + std::sqrt(std::pow(disk.b, 2) + z2), 2));
}

std::complex<double> Baryons_A::psi_MN_d2R2(std::complex<double> R2, std::complex<double> z2, const disk_3p& disk) {
    return disk.M / std::sqrt(R2 + std::pow(disk.a + std::sqrt(std::pow(disk.b, 2) + z2), 2));
}

std::complex<double> Baryons_A::psi_MN_d2z2(std::complex<double> R2, std::complex<double> z2, const disk_3p& disk) {
    return disk.M / std::sqrt(R2 + std::pow(disk.a + std::sqrt(std::pow(disk.b, 2) + z2), 2));
}
