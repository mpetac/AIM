#include <complex>
#include <iostream>

#include "src/Inversion.hpp"
#include "src/Model.hpp"
#include "src/Structs.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    
    halo_2p halo = {1e-2, 13};
    disk_3p disk1 = {1, 1, 0.1};
    disk_3p disk2 = {0.5, 2, 0.3};
    bulge_2p bulge = {1, 0.5};
    
    Halo_NFW DM(halo);
    Baryons_H_2MN baryons(disk1, disk2, bulge);
    
    Model model(&DM, &baryons);
    
    Inversion psdf(&model, 20, 5);
    
    return 0;
}
