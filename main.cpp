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
    
    halo_6p p_abc = {1e7, 13., 1., 3., 1.2, 1.};
    //Halo_gNFW halo(p_abc);
    Halo_sABC halo(p_abc);
    
    disk_3p disk1 = {5e10, 3.6, 0.5};
    disk_3p disk2 = {0., 2., 0.3};
    bulge_2p bulge = {1e11, 1.};
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    
    Model model(&halo, &baryons);
    
    Inversion psdf(&model, 100, 20);
    
    Observables obs(&model, &psdf);
    
    bool verbose = 0;
    int nPts = 100;
    double Rpts[nPts], zpts[nPts], result[nPts];
    double logRmin = -0.30103, logRmax = 3.;
    for (int i = 0; i < nPts; i++) {
        Rpts[i] = std::pow(10., logRmin + (logRmax - logRmin) * i / (nPts - 1));
        zpts[i] = 0;
        result[i] = 0;
    }
    obs.rho(nPts, Rpts, zpts, result);
    std::ofstream out_density("out/density.dat");
    for (int i = 0; i < nPts; i++) {
        double rho_true = std::real(model.rho(Rpts[i] * Rpts[i], 0, Rpts[i]));
        out_density << Rpts[i] << "\t" << result[i] << "\t" << rho_true << "\n";
        if (verbose) std::cout << "rho(" << Rpts[i] << "): " << result[i] << " / " << rho_true << " (" << result[i] / rho_true << ")" << std::endl;
    }
    out_density.close();
    
    
    int nVel = 100;
    double pv_mag[2 * nVel], pv_merid[2 * nVel], pv_azim[2 * nVel];
    
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
    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    
    /*
    */
    
    return 0;
}
