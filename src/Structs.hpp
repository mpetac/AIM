#pragma once

#include <complex>

/**
 * Struct for three-parameter component.
 * @param M Mass of the disk [M_sol]
 * @param a Scale length of the disk [kpc]
 * @param b Scale heigh of the disk [kpc]
 */

struct disk_3p {
    double M;
    double a;
    double b;
}; 

/**
 * Struct for two-parameter bulge.
 * @param M Mass of the bulge [M_sol]
 * @param a Scale parameter of the bulge [kpc]
 */

struct bulge_2p {
    double M;
    double a;
}; 

/**
 * Struct for two-parameter halo density profile.
 * @param rho_s Scale density of the halo [M_sol / kpc^3]
 * @param r_s Scale radius of the halo [kpc]
 */

struct halo_2p {
    double rho_s;
    double r_s;
};

/**
 * Struct for three-parameter halo density profile.
 * @param rho_s Scale density of the halo [M_sol / kpc^3]
 * @param r_s Scale radius of the halo [kpc]
 * @param gamma Central density slope of the halo
 */

struct halo_3p {
    double rho_s;
    double r_s;
    double gamma;
};

/**
 * Struct for four-parameter halo density profile.
 * @param rho_s Scale density of the halo [M_sol / kpc^3]
 * @param r_s Scale radius of the halo [kpc]
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
 * @param omega Angular velocity [Hz]
 * @param r_a Rotatinal scale radius [kpc]
 */

struct halo_rot_2p {
    double omega;
    double r_a;
};

/**
 * Auxilary struct for computing the inverse of gravitational potential
 */

struct psi_inverse_params {
    void *model;
    std::complex<double> xi;
    std::complex<double> R2;
};

/**
 * Auxilary struct for computing the inversion contour integrals
 */

struct inversion_params {
    void *model;
    void *z2interp;
    double E;
    double Lz;
    double Rc;
    double psiEnv;
    double h;
};
