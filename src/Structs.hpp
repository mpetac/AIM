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
 * Struct for general halo density profile.
 * @param rho_s Scale density of the halo [M_sol / kpc^3]
 * @param r_s Scale radius of the halo [kpc]
 * @param alpha Density slope parameter 1
 * @param beta Density slope parameter 2
 * @param gamma Density slope parameter 3
 * @param q2 Square of the flattening parameter
 */

struct halo_6p {
    double rho_s;
    double r_s;
    double alpha;
    double beta;
    double gamma;
    double q;
};

/**
 * Struct with halo rotation parameters
 * @param omega Angular velocity [Hz]
 * @param r_a Rotatinal scale radius [kpc]
 */

struct halo_rot_3p {
    double omega;
    double r_a;
    double q;
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

struct velocity_int_params {
    void *model;
    void *inversion;
    size_t nIntervals;
    double tolerance;
    double R;
    double psiRz;
    double v;
    double vf;
    int moment;
};

struct relative_velocity_int_params {
    void *model;
    void *inversion;
    size_t nIntervals;
    double tolerance;
    double R;
    double psiRz;
    double vMax;
    double v_rel;
    double uf;
    double um;
    double wf;
    int moment;
};

struct potential_int_params {
    void *halo;
    std::complex<double> R2;
    std::complex<double> z2;
    std::complex<double> e;
    std::complex<double> asin_e;
    double q;
    double tolerance;
    size_t nIntervals;
};

struct occupation_params {
    void *model;
    void *inversion;
    size_t nIntervals;
    double Emin;
    double Emax;
    double Lzmin;
    double Lzmax;
    double tolerance;
    double R;
    double z;
    double psiRz;
    double E;
};

struct dd_params {
    void *model;
    void *inversion;
    double cE;
    double sE;
    double power;
    double psiR;
    double R;
    double vR_sol;
    double vPhi_sol;
    double vz_sol;
    double vEarth;
    double vmin;
    double ct;
    double st;
    double cf;
    double sf;
    double tolerance;
    size_t nIntervals;
};
