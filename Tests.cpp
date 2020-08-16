#include <CppUTest/CommandLineTestRunner.h>
#include <CppUTest/TestHarness.h>
#include <CppUTest/PlatformSpecificFunctions.h>
#include <complex>

#include "Structs.hpp"
#include "Parametric_funcs.hpp"

TEST_GROUP(ParametricFuncsTestGroup) {
};

TEST(ParametricFuncsTestGroup, Potential_MyiamotoNagai) {
    
    disk_3p disk = {3.1, 1.2, 0.2};
    std::complex<double> R2(2.3, 1.7);
    std::complex<double> z2(0.4, 0.8);
    
    std::complex<double> val_psi = Parametric_funcs::psi_MN(R2, z2, disk);
    DOUBLES_EQUAL(4.80163e-6, std::real(val_psi), 1e-5);
    DOUBLES_EQUAL(-1.32217e-6, std::imag(val_psi), 1e-5);
    
    std::complex<double> val_psi_dR2 = Parametric_funcs::psi_MN_dR2(R2, z2, disk);
    DOUBLES_EQUAL(-2.40653e-7, std::real(val_psi_dR2), 1e-5);
    DOUBLES_EQUAL(2.50829e-7, std::imag(val_psi_dR2), 1e-5);
    
    std::complex<double> val_psi_dz2 = Parametric_funcs::psi_MN_dz2(R2, z2, disk);
    DOUBLES_EQUAL(-3.40482e-7, std::real(val_psi_dz2), 1e-5);
    DOUBLES_EQUAL(6.75805e-7, std::imag(val_psi_dz2), 1e-5);
    
    std::complex<double> val_psi_d2z2 = Parametric_funcs::psi_MN_d2z2(R2, z2, disk);
    DOUBLES_EQUAL(-3.01896e-7, std::real(val_psi_d2z2), 1e-5);
    DOUBLES_EQUAL(-4.81819e-7, std::imag(val_psi_d2z2), 1e-5);
    
    std::complex<double> val_psi_d2R2z2 = Parametric_funcs::psi_MN_d2R2z2(R2, z2, disk);
    DOUBLES_EQUAL(-1.11935e-8, std::real(val_psi_d2R2z2), 1e-5);
    DOUBLES_EQUAL(-1.58053e-7, std::imag(val_psi_d2R2z2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Potential_Hernquist) {
    
    bulge_2p bulge = {1.3, 0.4};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_psi = Parametric_funcs::psi_H(r, bulge);
    DOUBLES_EQUAL(4.80163e-6, std::real(val_psi), 1e-5);
    DOUBLES_EQUAL(-1.32217e-6, std::imag(val_psi), 1e-5);
    
    std::complex<double> val_psi_dr2 = Parametric_funcs::psi_H_dr2(r, bulge);
    DOUBLES_EQUAL(-2.40653e-7, std::real(val_psi_dr2), 1e-5);
    DOUBLES_EQUAL(2.50829e-7, std::imag(val_psi_dr2), 1e-5);
    
    std::complex<double> val_psi_d2r2 = Parametric_funcs::psi_H_d2r2(r, bulge);
    DOUBLES_EQUAL(-3.40482e-7, std::real(val_psi_d2r2), 1e-5);
    DOUBLES_EQUAL(6.75805e-7, std::imag(val_psi_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Potential_NFW) {
    
    halo_2p halo = {1.3, 0.4};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_psi = Parametric_funcs::psi_NFW(r, halo);
    DOUBLES_EQUAL(4.80163e-6, std::real(val_psi), 1e-5);
    DOUBLES_EQUAL(-1.32217e-6, std::imag(val_psi), 1e-5);
    
    std::complex<double> val_psi_dr2 = Parametric_funcs::psi_NFW_dr2(r, halo);
    DOUBLES_EQUAL(-2.40653e-7, std::real(val_psi_dr2), 1e-5);
    DOUBLES_EQUAL(2.50829e-7, std::imag(val_psi_dr2), 1e-5);
    
    std::complex<double> val_psi_d2r2 = Parametric_funcs::psi_NFW_d2r2(r, halo);
    DOUBLES_EQUAL(-3.40482e-7, std::real(val_psi_d2r2), 1e-5);
    DOUBLES_EQUAL(6.75805e-7, std::imag(val_psi_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Potential_BUR) {
    
    halo_2p halo = {1.3, 0.4};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_psi = Parametric_funcs::psi_BUR(r, halo);
    DOUBLES_EQUAL(4.80163e-6, std::real(val_psi), 1e-5);
    DOUBLES_EQUAL(-1.32217e-6, std::imag(val_psi), 1e-5);
    
    std::complex<double> val_psi_dr2 = Parametric_funcs::psi_BUR_dr2(r, halo);
    DOUBLES_EQUAL(-2.40653e-7, std::real(val_psi_dr2), 1e-5);
    DOUBLES_EQUAL(2.50829e-7, std::imag(val_psi_dr2), 1e-5);
    
    std::complex<double> val_psi_d2r2 = Parametric_funcs::psi_BUR_d2r2(r, halo);
    DOUBLES_EQUAL(-3.40482e-7, std::real(val_psi_d2r2), 1e-5);
    DOUBLES_EQUAL(6.75805e-7, std::imag(val_psi_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Rho_NFW) {
    
    halo_2p halo = {1.3, 0.4};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_rho = Parametric_funcs::rho_NFW(r, halo);
    DOUBLES_EQUAL(4.80163e-6, std::real(val_rho), 1e-5);
    DOUBLES_EQUAL(-1.32217e-6, std::imag(val_rho), 1e-5);
    
    std::complex<double> val_rho_dr2 = Parametric_funcs::rho_NFW_dr2(r, halo);
    DOUBLES_EQUAL(-2.40653e-7, std::real(val_rho_dr2), 1e-5);
    DOUBLES_EQUAL(2.50829e-7, std::imag(val_rho_dr2), 1e-5);
    
    std::complex<double> val_rho_d2r2 = Parametric_funcs::rho_NFW_d2r2(r, halo);
    DOUBLES_EQUAL(-3.40482e-7, std::real(val_rho_d2r2), 1e-5);
    DOUBLES_EQUAL(6.75805e-7, std::imag(val_rho_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Rho_BUR) {
    
    halo_2p halo = {1.3, 0.4};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_rho = Parametric_funcs::rho_BUR(r, halo);
    DOUBLES_EQUAL(4.80163e-6, std::real(val_rho), 1e-5);
    DOUBLES_EQUAL(-1.32217e-6, std::imag(val_rho), 1e-5);
    
    std::complex<double> val_rho_dr2 = Parametric_funcs::rho_BUR_dr2(r, halo);
    DOUBLES_EQUAL(-2.40653e-7, std::real(val_rho_dr2), 1e-5);
    DOUBLES_EQUAL(2.50829e-7, std::imag(val_rho_dr2), 1e-5);
    
    std::complex<double> val_rho_d2r2 = Parametric_funcs::rho_BUR_d2r2(r, halo);
    DOUBLES_EQUAL(-3.40482e-7, std::real(val_rho_d2r2), 1e-5);
    DOUBLES_EQUAL(6.75805e-7, std::imag(val_rho_d2r2), 1e-5);
    
};

int main(int ac, char** av) {
   return CommandLineTestRunner::RunAllTests(ac, av);
}
