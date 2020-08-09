#pragma once

struct disk_3p {
    double M;
    double a;
    double b;
}; 

struct bulge_2p {
    double M;
    double a;
    double b;
}; 

struct halo_2p {
    double rho_s;
    double r_s;
};

struct halo_3p {
    double rho_s;
    double r_s;
    double gamma;
};

struct halo_4p {
    double rho_s;
    double r_s;
    double gamma;
    double q;
};

