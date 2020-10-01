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
    time_t tStart = time(NULL);
    bool verbose = 0;
    
    /*
    double *params = new double[16];
    std::string line;
    std::ifstream paramsf("params.dat");
    if (paramsf.is_open()) {
        for (int i = 0; i < 16; i++) {
            getline(paramsf,line);
            params[i] = std::atof(line.c_str());
            if (verbose) std::cout << i << ", " << params[i] << ", " << line << '\n';
        }
        paramsf.close();
    }
    
    disk_3p disk1 = {params[0], params[1], params[2]};
    disk_3p disk2 = {params[3], params[4], params[5]};
    bulge_2p bulge = {params[6], params[7]};
    halo_2p p_nfw = {std::pow(10., params[8] + 9.), params[9]};
    halo_rot_2p rot_nfw = {params[13], params[14]};
    */
    
    disk_3p disk1 = {0., 3.6, 0.5};
    disk_3p disk2 = {0, 1., 1.};
    bulge_2p bulge = {0., 0.5};
    halo_2p p_nfw = {1e7, 13.};
    halo_rot_2p rot_nfw = {0, 1.};
    
    Halo_NFW halo(p_nfw, rot_nfw);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model model(&halo, &baryons);
    
    Inversion psdf(&model, 100, 10, 1e-3, 1);
    
    Observables obs(&model, &psdf, 1);
    
    /*
    int nPts = 100;
    double Rpts[nPts], zpts[nPts], result[nPts];
    double logRmin = -1, logRmax = 3.;
    for (int i = 0; i < nPts; i++) {
        Rpts[i] = std::pow(10., logRmin + (logRmax - logRmin) * i / (nPts - 1));
        zpts[i] = 0;
    }
    obs.rho(nPts, Rpts, zpts, result);
    std::ofstream out_density("out/rho.dat");
    for (int i = 0; i < nPts; i++) {
        double rho_true = std::real(model.rho(Rpts[i] * Rpts[i], 0, Rpts[i]));
        out_density << Rpts[i] << "\t" << result[i] << "\t" << rho_true << "\n";
        if (verbose) std::cout << "rho(" << Rpts[i] << "): " << result[i] << " / " << rho_true << " (" << result[i] / rho_true << ")" << std::endl;
    }
    out_density.close();
    
    std::ofstream out_v_mom("out/moments.dat");
    for (int i = 0; i < nPts; i++) {
        double mom_m2 = obs.v_mom(-2, Rpts[i], 0);
        double mom_m1 = obs.v_mom(-1, Rpts[i], 0);
        double mom_p1 = obs.v_mom(1, Rpts[i], 0);
        double mom_p2 = obs.v_mom(2, Rpts[i], 0);
        
        out_v_mom << Rpts[i] << "\t" << mom_m2 << "\t" << mom_m1 << "\t" << mom_p1 << "\t" << mom_p2 << "\n";
        if (verbose) std::cout << "v_mom:" << mom_m2 << ", " << mom_m1 << ", " << mom_p1 << ", " << mom_p2 << std::endl;
    }
    out_v_mom.close();
    */
    
    int nVel = 100;
    double pv_mag[2 * nVel], pv_merid[2 * nVel], pv_azim[2 * nVel], pv_rad[2 * nVel];
    
    obs.pv_mag(nVel, 8.122, 0, pv_mag);
    std::ofstream out_pv_mag("out/velocity_mag.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_mag << pv_mag[2 * i] << "\t" << pv_mag[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_mag(" << pv_mag[2 * i] << "): " << pv_mag[2 * i + 1] << std::endl;
    }
    out_pv_mag.close();
    
    
    obs.pv_merid(nVel, 8.122, 0, pv_merid);
    std::ofstream out_pv_merid("out/velocity_merid.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_merid << pv_merid[2 * i] << "\t" << pv_merid[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_merid(" << pv_merid[2 * i] << "): " << pv_merid[2 * i + 1] << std::endl;
    }
    out_pv_merid.close();
    
    obs.pv_azim(nVel, 8.122, 0, pv_azim);
    std::ofstream out_pv_azim("out/velocity_f.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_azim << pv_azim[2 * i] << "\t" << pv_azim[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_azim(" << pv_azim[2 * i] << "): " << pv_azim[2 * i + 1] << std::endl;
    }
    out_pv_azim.close();
    
    obs.pv_rad(nVel, 8.122, 0, pv_rad);
    std::ofstream out_pv_rad("out/velocity_R.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_rad << pv_rad[2 * i] << "\t" << pv_rad[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_rad(" << pv_rad[2 * i] << "): " << pv_rad[2 * i + 1] << std::endl;
    }
    out_pv_rad.close();

    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    
    return 0;
}
