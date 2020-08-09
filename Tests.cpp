#include <CppUTest/CommandLineTestRunner.h>
#include <CppUTest/TestHarness.h>
#include <CppUTest/PlatformSpecificFunctions.h>
#include <complex>

#include "Structs.hpp"
#include "Parametric_funcs.hpp"

TEST_GROUP(ParametricFuncsTestGroup) {
};

TEST(ParametricFuncsTestGroup, Potential_MN) {
    
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
    
}

int main(int ac, char** av) {
   return CommandLineTestRunner::RunAllTests(ac, av);
}
