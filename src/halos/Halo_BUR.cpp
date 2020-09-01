#include "Halo_BUR.hpp"

std::complex<double> Halo_BUR::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_dr2(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_dr2(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_d2r2(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::psi_NFW_d2r2(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::rho_NFW(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::rho_NFW_dr2(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Parametric_funcs::rho_NFW_d2r2(r, Halo_BUR::halo);
}

std::complex<double> Halo_BUR::v_phi(std::complex<double> R2) {
    return Parametric_funcs::v_phi(R2, Halo_BUR::halo_rot);
}


bool Halo_BUR::is_rotating() {
    if (Halo_BUR::halo_rot.omega == 0) return 0;
    else return 1;
}
