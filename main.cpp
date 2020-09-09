#include <complex>
#include <iostream>
#include <time.h>

#include "src/Observables.hpp"
#include "src/Inversion.hpp"
#include "src/Model.hpp"
#include "src/Structs.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    
    time_t tStart = time(NULL);
    
    halo_2p halo = {1e-2, 13};
    disk_3p disk1 = {1, 1, 0.1};
    disk_3p disk2 = {0.5, 2, 0.3};
    bulge_2p bulge = {1, 0.5};
    
    Halo_NFW DM(halo);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model model(&DM, &baryons);
    
    Inversion psdf(&model, 100, 10);
    
    Observables obs(&model, &psdf);
    
    int nPts = 10;
    double Rpts[nPts], zpts[nPts], result[nPts];
    double Rmin = 1e-3, Rmax = 1e2;
    for (int i = 0; i < nPts; i++) {
        //Rpts[i] = (Rmax - Rmin) * std::exp(
    }
    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    
    return 0;
}
