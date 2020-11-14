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

void print_density(int N, Model model, Observables obs, bool verbose, char *suffix) {
    double Rpts[N], zpts[N], result[N];
    double logRmin = -1, logRmax = 3.;
    for (int i = 0; i < N; i++) {
        Rpts[i] = std::pow(10., logRmin + (logRmax - logRmin) * i / (N - 1));
        zpts[i] = 0;
    }
    obs.rho(N, Rpts, zpts, result);
    
    char name[128];
    sprintf(name, "out/density_%s.dat", suffix);
    std::ofstream out_density(name);
    for (int i = 0; i < N; i++) {
        double rho_true = std::real(model.rho(Rpts[i] * Rpts[i], 0, Rpts[i]));
        out_density << Rpts[i] << "\t" << result[i] << "\t" << rho_true << "\n";
        if (verbose) std::cout << "rho(" << Rpts[i] << "): " << result[i] << " / " << rho_true << " (" << result[i] / rho_true << ")" << std::endl;
    }
    out_density.close();
}

void print_pv(int type, int N, double R, double z, Observables obs, bool verbose, char *suffix) {
    double pv[2 * N];
    char type_name[6];
    
    switch(type) {
        case 1:
            obs.pv_mag(N, R, z, pv);
            sprintf(type_name, "mag");
            break;
        case 2:
            obs.pv_merid(N, R, z, pv);
            sprintf(type_name, "merid");
            break;
        case 3:
            obs.pv_azim(N, R, z, pv);
            sprintf(type_name, "azim");
            break;
        case 4:
            obs.pv_rad(N, R, z, pv);
            sprintf(type_name, "rad");
            break;
        case 5:
            obs.pv_rel(N, R, z, pv);
            sprintf(type_name, "rel");
            break;
    }
    
    char name[128];
    sprintf(name, "out/pv_%s_%s.dat", type_name, suffix);
    std::ofstream out_pv(name);
    for (int i = 0; i < N; i++) {
        out_pv << pv[2 * i] << "\t" << pv[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_" << type_name << "(" << pv[2 * i] << "): " << pv[2 * i + 1] << std::endl;
    }
    out_pv.close();
}

void print_occupation(int N_E, int N_Lz, Observables obs, bool verbose, char *suffix) {
    double Epts[N_E], Lzpts[N_Lz], occ[(N_E - 1) * (N_Lz - 1)];
    for (int i = 0; i < N_E; i++) Epts[i] = 0.89 * i / (N_E - 1.) + 0.1;
    for (int i = 0; i < N_Lz; i++) Lzpts[i] = 2. * i / (N_Lz - 1.) - 1.;
    obs.occupation(N_E, N_Lz, Epts, Lzpts, occ, 1e-1);
    
    char name[128];
    sprintf(name, "out/occupation_%s.dat", suffix);
    std::ofstream out_occupation(name);
    for (int i = 0; i < N_E - 1; i++) {
        for (int j = 0; j < N_Lz - 1; j++) {
            out_occupation << occ[i * (N_Lz - 1) + j] << "\t";
            if (verbose) std::cout << "occupation(" << i << ", " << j << "):\t" << occ[i * (N_Lz - 1) + j] << std::endl;
        }
        out_occupation << std::endl;
    }
    out_occupation.close();
}

int main(int argc, char **argv) {
    //std::cout << "Hello, world!" << std::endl;
    
    time_t tStart = time(NULL);
    
    // Define struct with DM halo parameters. Here we asume spherical DM density profile with density 1e7 M_sol / kpc^3 and scale density of 13 kpc.
    halo_2p p_nfw = {1e7, 13.};
    // Initialize the DM halo object.
    Halo_NFW halo(p_nfw);
    
    // Define structs related to the baryonic distribution. In this example we assume a model consisting of tow Myiamoto-Nagai disks and a spherical Hernquist bulge.
    disk_3p disk1 = {5e10, 3.6, 0.4};
    disk_3p disk2 = {0., 1., 1.};
    bulge_2p bulge = {7e9, 0.5};
    // Initialize the baryonic model
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    // Initialize the galactic model using the previously defined halo and baryons
    Model model(&halo, &baryons);
    
    // Interpolate the PSDF obtained for the specified galactic model with given number of relative energy and angular momentum points
    Inversion psdf(&model, 500, 20, 1e-3, 1);
    
    // Initialize the class for computing various observable quantities from the PSDF (namely DM density and various projections of the velocity distribution)
    Observables obs(&model, &psdf);
    
    
    // Tabulate the DM density distribution and write it to a file
    bool verbose = 0;
    double R = 8.122, z = 0;
    char suffix[64];
    sprintf(suffix, "test");
    
    std::cout << obs.v_rel_mom(0, 100., z, 1e-1) << std::endl;
    
    /*
    print_density(50, model, obs, verbose, suffix);
    print_pv(1, 100, R, z, obs, verbose, suffix);
    print_pv(2, 100, R, z, obs, verbose, suffix);
    print_pv(3, 100, R, z, obs, verbose, suffix);
    print_pv(4, 100, R, z, obs, verbose, suffix);
    //print_pv(5, 7, R, z, obs, verbose, suffix);
    print_occupation(10, 5, obs, verbose, suffix);
    */
    
    /* * * VELOCITY DISTRIBUTIONS * * */
    
    /*
    int nVel = 30;
    double pv_mag[2 * nVel], pv_merid[2 * nVel], pv_azim[2 * nVel], pv_rad[2 * nVel], pv_rel[2 * nVel];
    
    // Tabulate the magnitude of DM velocity distribution and write it to a file
    obs.pv_mag(nVel, 8.122, 0, pv_mag);
    std::ofstream out_pv_mag("out/pv_mag.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_mag << pv_mag[2 * i] << "\t" << pv_mag[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_mag(" << pv_mag[2 * i] << "): " << pv_mag[2 * i + 1] << std::endl;
    }
    out_pv_mag.close();
    
    // Tabulate the DM velocity distribution in meridional plane and write it to a file
    obs.pv_merid(nVel, 8.122, 0, pv_merid);
    std::ofstream out_pv_merid("out/pv_merid.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_merid << pv_merid[2 * i] << "\t" << pv_merid[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_merid(" << pv_merid[2 * i] << "): " << pv_merid[2 * i + 1] << std::endl;
    }
    out_pv_merid.close();
    
    // Tabulate the DM velocity distribution in azimuthal direction and write it to a file
    obs.pv_azim(nVel, 8.122, 0, pv_azim);
    std::ofstream out_pv_azim("out/pv_azim.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_azim << pv_azim[2 * i] << "\t" << pv_azim[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_azim(" << pv_azim[2 * i] << "): " << pv_azim[2 * i + 1] << std::endl;
    }
    out_pv_azim.close();
    
    // Tabulate the DM velocity distribution in radial direction and write it to a file
    obs.pv_rad(nVel, 8.122, 0, pv_rad);
    std::ofstream out_pv_rad("out/pv_rad.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_rad << pv_rad[2 * i] << "\t" << pv_rad[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_rad(" << pv_rad[2 * i] << "): " << pv_rad[2 * i + 1] << std::endl;
    }
    out_pv_rad.close();
    
    
    obs.pv_rel(nVel, 8.122, 0, pv_rel, 1e-1);
    std::ofstream out_pv_rel("out/pv_rel.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_rel << pv_rel[2 * i] << "\t" << pv_rel[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_rel(" << pv_rel[2 * i] << "): " << pv_rel[2 * i + 1] << std::endl;
    }
    out_pv_rel.close();
    */
    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    return 0;
}
