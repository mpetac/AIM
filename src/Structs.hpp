#pragma once

#include <complex>

/**
 * Struct for three-parameter component.
 * @param M Mass of the disk
 * @param a Scale length of the disk
 * @param b Scale heigh of the disk
 */

struct disk_3p {
    double M;
    double a;
    double b;
}; 

/**
 * Struct for two-parameter bulge.
 * @param M Mass of the bulge
 * @param a Scale parameter of the bulge
 */

struct bulge_2p {
    double M;
    double a;
}; 

/**
 * Struct for two-parameter halo density profile.
 * @param rho_s Scale density of the halo
 * @param r_s Scale radius of the halo
 */

struct halo_2p {
    double rho_s;
    double r_s;
};

/**
 * Struct for three-parameter halo density profile.
 * @param rho_s Scale density of the halo
 * @param r_s Scale radius of the halo
 * @param gamma Central density slope of the halo
 */

struct halo_3p {
    double rho_s;
    double r_s;
    double gamma;
};

/**
 * Struct for four-parameter halo density profile.
 * @param rho_s Scale density of the halo
 * @param r_s Scale radius of the halo
 * @param gamma Central density slope of the halo
 * @param q Flattening of the halo along the axis of symmetry
 */

struct halo_4p {
    double rho_s;
    double r_s;
    double gamma;
    double q;
};

/**
 * Struct with halo rotation parameters
 * @param omega Angular velocity
 * @param r_a Rotatinal scale radius
 */

struct halo_rot_2p {
    double omega;
    double r_a;
};

struct psi_inverse_params {
    void *model;
    std::complex<double> xi;
    std::complex<double> R2;
};

struct inversion_params {
    void *model;
    void *z2interp;
    double E;
    double Lz;
    double Rc;
    double psiEnv;
    double h;
};
