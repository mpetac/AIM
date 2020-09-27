#include <complex>
#include <iostream>
#include <fstream>
#include <time.h>

#include "src/Observables.hpp"
#include "src/Inversion.hpp"
#include "src/Model.hpp"
#include "src/Structs.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/halos/Halo_gNFW.hpp"
#include "src/halos/Halo_sABC.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    
    time_t tStart = time(NULL);
    
//     halo_2p p_nfw = {1e7, 13.};
    halo_2p p_nfw = {4.02e5, 49.2};
    Halo_NFW halo(p_nfw);
    
//     halo_6p p_abc = {4.02e5, 49.2, 1., 3., 1.01, 1.};
//     Halo_gNFW halo(p_abc);
//     Halo_sABC halo(p_abc);
    
    disk_3p disk1 = {1.39e10, 20., 10.};
    disk_3p disk2 = {4.12e10, 9.52, 0.};
    bulge_2p bulge = {1.59e10, 0.484};
//     disk_3p disk1 = {5e10, 3.6, 0.5};
//     disk_3p disk2 = {0., 9.52, 0.};
//     bulge_2p bulge = {1e11, 1.};
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model model(&halo, &baryons);
    
    /*
    int N = 100;
    std::ofstream out("out/d2rho_dpsi2.dat");
    for (int i = 0; i < N; i++) {
        std::complex<double> R2(1000., 0.);
        std::complex<double> z2 = 1000. * i / (N - 1.);
        std::complex<double> r = std::sqrt(R2 + z2);
        out << std::real(z2) << "\t" << std::real(model.rho_d2psi2(R2, z2, r)) << std::endl;
    }
    out.close();
    */
    
    Inversion psdf(&model, 1000, 10);
    
    Observables obs(&model, &psdf);
    
    
    
    bool verbose = 0;
    int nPts = 20;
    double Rpts[nPts], zpts[nPts], result[nPts];
    double logRmin = -1, logRmax = 3.;
    for (int i = 0; i < nPts; i++) {
        Rpts[i] = std::pow(10., logRmin + (logRmax - logRmin) * i / (nPts - 1));
        zpts[i] = 0;
    }
    obs.rho(nPts, Rpts, zpts, result);
    std::ofstream out_density("out/density.dat");
    for (int i = 0; i < nPts; i++) {
        double rho_true = std::real(model.rho(Rpts[i] * Rpts[i], 0, Rpts[i]));
        out_density << Rpts[i] << "\t" << result[i] << "\t" << rho_true << "\n";
        if (verbose) std::cout << "rho(" << Rpts[i] << "): " << result[i] << " / " << rho_true << " (" << result[i] / rho_true << ")" << std::endl;
    }
    out_density.close();
    
    std::ofstream out_v_mom("out/v_mom.dat");
    for (int i = 0; i < nPts; i++) {
        double mom_m2 = obs.v_mom(-2, Rpts[i], 0);
        double mom_m1 = obs.v_mom(-1, Rpts[i], 0);
        double mom_p1 = obs.v_mom(1, Rpts[i], 0);
        double mom_p2 = obs.v_mom(2, Rpts[i], 0);
        
        out_v_mom << Rpts[i] << "\t" << mom_m2 << "\t" << mom_m1 << "\t" << mom_p1 << "\t" << mom_p2 << "\n";
        if (verbose) std::cout << "v_mom:" << mom_m2 << ", " << mom_m1 << ", " << mom_p1 << ", " << mom_p2 << std::endl;
    }
    out_v_mom.close();
    
    
    int nVel = 100;
    double pv_mag[2 * nVel], pv_merid[2 * nVel], pv_azim[2 * nVel], pv_rad[2 * nVel];
    
    obs.pv_mag(nVel, 8.122, 0, pv_mag);
    std::ofstream out_pv_mag("out/pv_mag.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_mag << pv_mag[2 * i] << "\t" << pv_mag[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_mag(" << pv_mag[2 * i] << "): " << pv_mag[2 * i + 1] << std::endl;
    }
    out_pv_mag.close();
    
    obs.pv_merid(nVel, 8.122, 0, pv_merid);
    std::ofstream out_pv_merid("out/pv_merid.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_merid << pv_merid[2 * i] << "\t" << pv_merid[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_merid(" << pv_merid[2 * i] << "): " << pv_merid[2 * i + 1] << std::endl;
    }
    out_pv_merid.close();
    
    obs.pv_azim(nVel, 8.122, 0, pv_azim);
    std::ofstream out_pv_azim("out/pv_azim.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_azim << pv_azim[2 * i] << "\t" << pv_azim[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_azim(" << pv_azim[2 * i] << "): " << pv_azim[2 * i + 1] << std::endl;
    }
    out_pv_azim.close();
    
    obs.pv_rad(nVel, 8.122, 0, pv_rad);
    std::ofstream out_pv_rad("out/pv_rad.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_rad << pv_rad[2 * i] << "\t" << pv_rad[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_rad(" << pv_rad[2 * i] << "): " << pv_rad[2 * i + 1] << std::endl;
    }
    out_pv_rad.close();
    
    /* */
    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    
    return 0;
}
