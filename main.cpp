#include <iostream>
#include "Baryons_A.hpp"
#include "Structs.hpp"

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    
    disk_3p disk1 = {1, 1, 0.1};
    disk_3p disk2 = {0.5, 2, 0.3};
    bulge_2p bulge = {1, 0.5};
    
    Baryons_A bar(disk1, disk2, bulge);
    
    std::cout << bar.psi(1, 1) << std::endl;
    
    return 0;
}
