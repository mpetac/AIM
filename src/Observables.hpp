#pragma once

#include <time.h>
#include <complex>
#include <future>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "Inversion.hpp"
#include "Model.hpp"
#include "Structs.hpp"

/**
 * Class used for computing observables from the phase-space distribution function (PSDF).
 */ 

class Observables {
    
private:
    /// Galactic model
    Model *model;
    /// Instance of the inversion class
    Inversion *inversion;
    /// Number of intervals used in the numerical integration
    size_t nIntervals = 1e5;
    /// Parameter controling verbose output
    bool verbose;
    /// Maximum velocity to which the velocity probability distribution is computed
    double vMax = 1000.;
    
    /// Function for computing the DM density at a given point
    double rho_int(double R, double z, double tolerance);
    
    /// Function for computing the DM's meridional velocity distribution at a given point
    double pv_merid_int(double v_merid, double R, double psiRz, double tolerance);
    
    /// Function for computing the DM's azimuthal velocity distribution at a given point
    double pv_azim_int(double v_merid, double R, double psiRz, double tolerance);
    
    /// Function for computing the DM's radial (or equivalently along z-direction) velocity distribution at a given point
    double pv_rad_int(double v_merid, double R, double psiRz, double tolerance);
    
    /// Function for computing the DM's velocity magnitude distribution at a given point
    double pv_rel_int(double v_rel, double R, double psiRz, double tolerance);
    
    /// Function for computing the mass in a given phase-space region
    double occupation_int(double Emin, double Emax, double Lzmin, double Lzmax, double tolerance);
    
    /// Function for computing the DD astrophysical factor
    double dd_int(double vmin, dd_params params);
    
    /// GSL error handler
    static void GSL_error_func(const char * reason, const char * file, int line, int gsl_errno);
    /// Silent GSL error handler
    static void GSL_error_func_silent(const char * reason, const char * file, int line, int gsl_errno);
    
public:
    /// Initializer
    Observables(Model *model, Inversion *inversion, bool verbose = 0);

    /// Computes the DM density profile
    void rho(int N, double *Rpts, double *zpts, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity magnitude
    void pv_mag(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity in meridional plane
    void pv_merid(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity along the azimuthal direction
    void pv_azim(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity along the radial (or equivalently z) direction
    void pv_rad(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM relative velocity magnitude
    void pv_rel(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes moments of the DM's velocity distribution
    double v_mom(int mom, double R, double z, double tolerance=1e-3);
    
    /// Computes moments of the DM's relative velocity distribution
    double v_rel_mom(int mom, double R, double z, double tolerance=1e-3);
    
    /// Computes the mass in phase-space bins
    void occupation(int N_E, int N_Lz, double *Epts, double *Lzpts, double *result, double tolerance=1e-3);
    
    /// Computes the direct detection astrophysical factor
    void dd(int N, double *result, double t=0, double power=-1., double R=8.122, double vR=11., double vPhi=-1., double vz=7., double vEarth=30., double vmax=1000., double tolerance=1e-3);
    
};
