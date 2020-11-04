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

void print_file(char *name, double *results, int N, int M) {
    std::ofstream out_density(name);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            out_density << results[i * M + j] << "\t";
        }
        out_density << std::endl;
    }
    out_density.close();
}

int main(int argc, char **argv) {
    //std::cout << "Hello, world!" << std::endl;
    
    time_t tStart = time(NULL);
    
    // Define struct with DM halo parameters. Here we asume spherical DM density profile with density 1e7 M_sol / kpc^3 and scale density of 13 kpc.
    halo_2p p_nfw = {1e8, 13.};
    // Initialize the DM halo object.
    Halo_NFW halo(p_nfw);
    
    // Define structs related to the baryonic distribution. In this example we assume a model consisting of tow Myiamoto-Nagai disks and a spherical Hernquist bulge.
    disk_3p disk1 = {0, 3.6, 0.4};
    disk_3p disk2 = {0., 1., 1.};
    bulge_2p bulge = {7e9, 0.5};
    // Initialize the baryonic model
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    // Initialize the galactic model using the previously defined halo and baryons
    Model model(&halo, &baryons);
    
    // Interpolate the PSDF obtained for the specified galactic model with given number of relative energy and angular momentum points
    Inversion psdf(&model, 100, 10, 1e-3, 1);
    
    // Initialize the class for computing various observable quantities from the PSDF (namely DM density and various projections of the velocity distribution)
    Observables obs(&model, &psdf);
    
    //double R = 250.;
    //std::cout << -2. * std::pow(R, 3) * std::real(model.psi_dR2(R * R, 0, R)) / 4.3e-6 << ", " << std::real(model.psi(R * R, 0, R)) / model.psi0 << std::endl;
    
    double Mtot = obs.occupation_int(0.1, 1., -1., 1., 1e-3);
    std::cout << Mtot << std::endl;
    
    
    // Tabulate the DM density distribution and write it to a file
    bool verbose = 0;
    
    
    /* * * DENSITY * * */
    
    /*
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
    */
    
    /* * * OCCUPATION NUMBERS * * */
    
    /*
    int N_E = 10, N_Lz = 10;
    int Nocc = N_E * N_Lz;
    double occ[Nocc];
    obs.occupation(N_E, N_Lz, occ);
    std::ofstream out_occupation("out/occupation.dat");
    for (int i = 0; i < N_E; i++) {
        for (int j = 0; j < N_Lz; j++) {
            out_occupation << occ[i * N_Lz + j] << "\t";
            if (verbose) std::cout << "occupation(" << i << ", " << j << "):\t" << occ[i * N_Lz + j] << std::endl;
        }
        out_occupation << std::endl;
    }
    out_occupation.close();
    */
    
    /* * * VELOCITY DISTRIBUTIONS * * */
    
    /*
    int nVel = 100;
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
    */
    
    /*
    obs.pv_rel(nVel, 8.122, 0, pv_rel);
    std::ofstream out_pv_rel("out/pv_rel.dat");
    for (int i = 0; i < nVel; i++) {
        out_pv_rel << pv_rel[2 * i] << "\t" << pv_rel[2 * i + 1] << "\n";
        if (verbose) std::cout << "pv_rel(" << pv_rel[2 * i] << "): " << pv_rel[2 * i + 1] << std::endl;
    }
    out_pv_rel.close();
    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    */
    return 0;
}
