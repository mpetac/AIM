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
    
    halo_2p halo = {1e7, 13.};
    disk_3p disk1 = {0., 3.6, 0.5};
    disk_3p disk2 = {0., 2., 0.3};
    bulge_2p bulge = {0., 1.};
    
    Halo_NFW DM(halo);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model model(&DM, &baryons);
    
    
    Inversion psdf(&model, 100, 20);
    
    //psdf.test();
    
    /*
    int N = 100;
    for (int i = 0; i < N; i++) {
        double E = 1. - pow(10., -4. + 4. * i / (N - 1.));
        double F = psdf.eval_F(E, 0);
        std::cout << "F(" << E << ", 0) = " << F << std::endl;
    }
    */
    
    
    Observables obs(&model, &psdf);
    
    int nPts = 5;
    double Rpts[nPts], zpts[nPts], result[nPts];
    double logRmin = -0.30103, logRmax = 3.;
    for (int i = 0; i < nPts; i++) {
        Rpts[i] = std::pow(10., logRmin + (logRmax - logRmin) * i / (nPts - 1));
        zpts[i] = 0;
        result[i] = 0;
    }
    
    obs.rho(nPts, Rpts, zpts, result);
    
    for (int i = 0; i < nPts; i++) {
        double rho_true = std::real(model.rho(Rpts[i] * Rpts[i], 0, Rpts[i]));
        std::cout << "rho(" << Rpts[i] << "): " << result[i] << " / " << rho_true << " (" << result[i] / rho_true << ")" << std::endl;
    }
    
    
    double dt = difftime(time(NULL), tStart);
    std::cout << "Done in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    
    
    return 0;
}
