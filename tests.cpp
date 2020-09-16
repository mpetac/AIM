#include <CppUTest/CommandLineTestRunner.h>
#include <CppUTest/TestHarness.h>
#include <CppUTest/PlatformSpecificFunctions.h>
#include <complex>
#include <iostream>

#include "src/Structs.hpp"
#include "src/Parametric_funcs.hpp"
#include "src/Numeric_funcs.hpp"
#include "src/Model.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/halos/Halo_gNFW.hpp"
#include "src/halos/Halo_sABC.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"
class Halo_sABC;

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

TEST(ParametricFuncsTestGroup, Rho_sABC) {
    halo_6p halo = {1e-2, 16., 1.8, 3.1, 1.2, 0.9};
    std::complex<double> R2(1.3, 2.7);
    std::complex<double> z2(2.8, 0.4);
    std::complex<double> m2 = R2 + z2 / std::pow(halo.q, 2);
    /*
    halo_2p halo_nfw = {1e-2, 16.};
    std::complex<double> r = std::sqrt(R2 + z2);
    std::cout << Parametric_funcs::rho_NFW(r, halo_nfw) << std::endl;*/
    
    std::complex<double> val_rho = Parametric_funcs::rho_sABC(m2, halo);
    DOUBLES_EQUAL(0.0884148, std::real(val_rho), 1e-5);
    DOUBLES_EQUAL(-0.034478, std::imag(val_rho), 1e-5);
    
    std::complex<double> val_rho_dr2 = Parametric_funcs::rho_sABC_dz2(m2, halo);
    DOUBLES_EQUAL(-0.00756271, std::real(val_rho_dr2), 1e-5);
    DOUBLES_EQUAL(0.0103398, std::imag(val_rho_dr2), 1e-5);
    
    std::complex<double> val_rho_d2r2 = Parametric_funcs::rho_sABC_d2z2(m2, halo);
    DOUBLES_EQUAL(0.000164061, std::real(val_rho_d2r2), 1e-5);
    DOUBLES_EQUAL(-0.00438241, std::imag(val_rho_d2r2), 1e-5);
    
};

TEST(ParametricFuncsTestGroup, d2Rho_dPsi2) {
    halo_2p halo = {1e7, 13};
    disk_3p disk1 = {5e10, 3.6, 0.5};
    disk_3p disk2 = {0., 2., 0.3};
    bulge_2p bulge = {1e11, 1.};
    
    Halo_NFW DM(halo);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    Model model(&DM, &baryons);
    
    std::complex<double> R2(7.3, -2.7);
    std::complex<double> z2(100.3, 5.7);
    std::complex<double> d2rho_dpsi2 = model.rho_d2psi2(R2, z2, std::sqrt(R2 + z2));
    DOUBLES_EQUAL(0.00218314, std::real(d2rho_dpsi2), 1e-5);
    DOUBLES_EQUAL(1.31776e-5, std::imag(d2rho_dpsi2), 1e-5);
};

TEST_GROUP(NumericFuncsTestGroup) {
};

TEST(NumericFuncsTestGroup, PsiSpheroid) {
    halo_2p p_nfw = {1e7, 13.};
    Halo_NFW nfw(p_nfw);
    
    halo_6p p_abc = {1e7, 13., 1., 3., 1., 1.};
    Halo_sABC abc(p_abc);
    
    double psi0_ana = std::real(nfw.psi(0., 0., 0.));
    double psi0_num = std::real(abc.psi(0., 0., 0.));
    DOUBLES_EQUAL(1., psi0_num / psi0_ana, 1e-3);
    
    std::complex<double> R2(7.3, -2.7);
    std::complex<double> z2(1.3, 5.7);
    std::complex<double> r = std::sqrt(R2 + z2);
    
    std::complex<double> psi_ana = nfw.psi(R2, z2, r);
    std::complex<double> psi_num = abc.psi(R2, z2, r);
    DOUBLES_EQUAL(1., std::real(psi_num / psi_ana), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num / psi_ana), 1e-3);
    
    std::complex<double> psi_ana_dR2 = nfw.psi_dR2(R2, z2, r);
    std::complex<double> psi_num_dR2 = abc.psi_dR2(R2, z2, r);
    DOUBLES_EQUAL(1., std::real(psi_num_dR2 / psi_ana_dR2), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num_dR2 / psi_ana_dR2), 1e-3);
    
    std::complex<double> psi_ana_dz2 = nfw.psi_dz2(R2, z2, r);
    std::complex<double> psi_num_dz2 = abc.psi_dz2(R2, z2, r);
    DOUBLES_EQUAL(1., std::real(psi_num_dz2 / psi_ana_dz2), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num_dz2 / psi_ana_dz2), 1e-3);
    
    std::complex<double> psi_ana_d2z2 = nfw.psi_d2z2(R2, z2, r);
    std::complex<double> psi_num_d2z2 = abc.psi_d2z2(R2, z2, r);
    DOUBLES_EQUAL(1., std::real(psi_num_d2z2 / psi_ana_d2z2), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num_d2z2 / psi_ana_d2z2), 1e-3);
    
    std::complex<double> psi_ana_d2R2z2 = nfw.psi_d2R2z2(R2, z2, r);
    std::complex<double> psi_num_d2R2z2 = abc.psi_d2R2z2(R2, z2, r);
    DOUBLES_EQUAL(1., std::real(psi_num_d2R2z2 / psi_ana_d2R2z2), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num_d2R2z2 / psi_ana_d2R2z2), 1e-3);
}

TEST(NumericFuncsTestGroup, PsigNFW) {
    halo_6p p_sabc = {1e7, 13., 1., 3., 1.01, 1.};
    Halo_sABC sabc(p_sabc);
    
    halo_6p p_gnfw = {1e7, 13., 1., 3., 1.01, 1.};
    Halo_gNFW gnfw(p_gnfw);
    
    std::complex<double> R2(17.3, -2.7);
    std::complex<double> z2(1.3, 10.7);
    std::complex<double> r = std::sqrt(R2 + z2);
    
    std::complex<double> psi_ana = sabc.psi(R2, z2, r);
    std::complex<double> psi_num = gnfw.psi(R2, z2, r);
    std::cout << psi_ana << ", " << psi_num << std::endl;
    DOUBLES_EQUAL(1., std::real(psi_num / psi_ana), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num / psi_ana), 1e-3);
    
    std::complex<double> psi_ana_dz2 = sabc.psi_dz2(R2, z2, r);
    std::complex<double> psi_num_dz2 = gnfw.psi_dz2(R2, z2, r);
    std::cout << psi_ana_dz2 << ", " << psi_num_dz2 << std::endl;
    DOUBLES_EQUAL(1., std::real(psi_num_dz2 / psi_ana_dz2), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num_dz2 / psi_ana_dz2), 1e-3);
    
    std::complex<double> psi_ana_d2z2 = sabc.psi_d2z2(R2, z2, r);
    std::complex<double> psi_num_d2z2 = gnfw.psi_d2z2(R2, z2, r);
    std::cout << psi_ana_d2z2 << ", " << psi_num_d2z2 << std::endl;
    DOUBLES_EQUAL(1., std::real(psi_num_d2z2 / psi_ana_d2z2), 1e-3);
    DOUBLES_EQUAL(0., std::imag(psi_num_d2z2 / psi_ana_d2z2), 1e-3);
};

TEST_GROUP(ModelTestGroup) {
};

TEST(ModelTestGroup, Rcirc) {
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
};

TEST(ModelTestGroup, PsiInverse) {
    halo_2p halo = {1e7, 13.};
    disk_3p disk1 = {5e10, 3.6, 0.5};
    disk_3p disk2 = {0, 1., 1.};
    bulge_2p bulge = {1e11, 1.};
    
    Halo_NFW DM(halo);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model m(&DM, &baryons);
    
    double E = 0.99 * std::real(m.psi(0, 0, 0));
    double Rc_E = m.Rcirc(E);
    double Lz = 1. * std::pow(Rc_E, 2) * std::sqrt(-2 * std::real(m.psi_dR2(std::pow(Rc_E, 2), 0, Rc_E)));
    std::complex<double> R2(1.18, 1.3);
    std::complex<double> xi = std::pow(Lz, 2) / (2. * R2) + E;
    std::complex<double> z2 = m.psi_inverse(xi, E, Lz);
    double rel_dif = std::abs(xi - m.psi(R2, z2, std::sqrt(R2 + z2))) / m.psi0;
    //std::cout << xi << " : " << m.psi(R2, z2, std::sqrt(R2 + z2)) << std::endl;
    DOUBLES_EQUAL(0., rel_dif, 1e-5);
};

int main(int ac, char** av) {
   return CommandLineTestRunner::RunAllTests(ac, av);
}
