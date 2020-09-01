#include <CppUTest/CommandLineTestRunner.h>
#include <CppUTest/TestHarness.h>
#include <CppUTest/PlatformSpecificFunctions.h>
#include <complex>
#include <iostream>

#include "src/Structs.hpp"
#include "src/Parametric_funcs.hpp"
#include "src/Model.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"

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
    
    halo_2p halo = {1e-2, 16};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_psi = Parametric_funcs::psi_NFW(r, halo);
    DOUBLES_EQUAL(1.31913e-4, std::real(val_psi), 1e-5);
    DOUBLES_EQUAL(-1.03862e-5, std::imag(val_psi), 1e-5);
    
    
    std::complex<double> val_psi_dr2 = Parametric_funcs::psi_NFW_dr2(r, halo);
    DOUBLES_EQUAL(-1.52629e-7, std::real(val_psi_dr2), 1e-5);
    DOUBLES_EQUAL(6.2169e-7, std::imag(val_psi_dr2), 1e-5);
    
    std::complex<double> val_psi_d2r2 = Parametric_funcs::psi_NFW_d2r2(r, halo);
    DOUBLES_EQUAL(-3.93567e-8, std::real(val_psi_d2r2), 1e-5);
    DOUBLES_EQUAL(1.06767e-8, std::imag(val_psi_d2r2), 1e-5);
};

TEST(ParametricFuncsTestGroup, Potential_BUR) {
    
    halo_2p halo = {1e-2, 16};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_psi = Parametric_funcs::psi_BUR(r, halo);
    DOUBLES_EQUAL(1.09075e-4, std::real(val_psi), 1e-5);
    DOUBLES_EQUAL(-6.49125e-7, std::imag(val_psi), 1e-5);
    
    std::complex<double> val_psi_dr2 = Parametric_funcs::psi_BUR_dr2(r, halo);
    DOUBLES_EQUAL(-8.45544e-8, std::real(val_psi_dr2), 1e-5);
    DOUBLES_EQUAL(1.14395e-8, std::imag(val_psi_dr2), 1e-5);
    
    std::complex<double> val_psi_d2r2 = Parametric_funcs::psi_BUR_d2r2(r, halo);
    DOUBLES_EQUAL(3.10035e-10, std::real(val_psi_d2r2), 1e-5);
    DOUBLES_EQUAL(-6.43393e-10, std::imag(val_psi_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Rho_NFW) {
    
    halo_2p halo = {1e-2, 16.};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_rho = Parametric_funcs::rho_NFW(r, halo);
    DOUBLES_EQUAL(0.00618087, std::real(val_rho), 1e-5);
    DOUBLES_EQUAL(-0.0441534, std::imag(val_rho), 1e-5);
    
    std::complex<double> val_rho_dr2 = Parametric_funcs::rho_NFW_dr2(r, halo);
    DOUBLES_EQUAL(0.00290775, std::real(val_rho_dr2), 1e-5);
    DOUBLES_EQUAL(-0.000907542, std::imag(val_rho_dr2), 1e-5);
    
    std::complex<double> val_rho_d2r2 = Parametric_funcs::rho_NFW_d2r2(r, halo);
    DOUBLES_EQUAL(0.000405253, std::real(val_rho_d2r2), 1e-5);
    DOUBLES_EQUAL(0.000303309, std::imag(val_rho_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, Rho_BUR) {
    
    halo_2p halo = {1e-2, 16.};
    std::complex<double> r(1.3, 2.7);
    
    std::complex<double> val_rho = Parametric_funcs::rho_BUR(r, halo);
    DOUBLES_EQUAL(0.00918295, std::real(val_rho), 1e-5);
    DOUBLES_EQUAL(-0.00169805, std::imag(val_rho), 1e-5);
    
    std::complex<double> val_rho_dr2 = Parametric_funcs::rho_BUR_dr2(r, halo);
    DOUBLES_EQUAL(-0.0000463169, std::real(val_rho_dr2), 1e-5);
    DOUBLES_EQUAL(0.0000962455, std::imag(val_rho_dr2), 1e-5);
    
    std::complex<double> val_rho_d2r2 = Parametric_funcs::rho_BUR_d2r2(r, halo);
    DOUBLES_EQUAL(-5.40276e-6, std::real(val_rho_d2r2), 1e-5);
    DOUBLES_EQUAL(1.19492e-6, std::imag(val_rho_d2r2), 1e-5);
    
};

TEST_GROUP(ModelTestGroup) {
};

TEST(ModelTestGroup, Model) {
    halo_2p halo = {1e-2, 13};
    disk_3p disk1 = {1, 1, 0.1};
    disk_3p disk2 = {0.5, 2, 0.3};
    bulge_2p bulge = {1, 0.5};
    
    Halo_NFW DM(halo);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model m(&DM, &baryons);
    
    double Rc = 2.13;
    double E_Rc = std::real(m.psi(std::pow(Rc, 2), 0, Rc) + std::pow(Rc, 2) * m.psi_dR2(std::pow(Rc, 2), 0, Rc));
    DOUBLES_EQUAL(m.Rcirc(E_Rc), Rc, 1e-3);
    
    
    double E = 0.78 * std::real(m.psi(0, 0, 0));
    double Rc_E = m.Rcirc(E);
    double Lz = 0.28 * std::real(m.psi_dR2(std::pow(Rc_E, 2), 0, Rc_E));
    std::complex<double> R2(1.18, 0.3);
    std::complex<double> xi = std::pow(Lz, 2) / (2. * R2) + E;
    std::complex<double> z2 = m.psi_inverse(xi, E, Lz);
    DOUBLES_EQUAL(std::real(xi), std::real(m.psi(R2, z2, std::sqrt(R2 + z2))), 1e-5);
    DOUBLES_EQUAL(std::imag(xi), std::imag(m.psi(R2, z2, std::sqrt(R2 + z2))), 1e-5);
};

int main(int ac, char** av) {
   return CommandLineTestRunner::RunAllTests(ac, av);
}
