#pragma once

#include <complex>
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
    
    /// Function for computing the DM density at a given point
    double rho_int(double R, double z, double tolerance);
    
public:
    /// Initializer
    Observables(Model *model, Inversion *inversion);
    
    /// Destructor
    ~Observables();

    /// Computes the DM density profile
    void rho(int N, double *Rpts, double *zpts, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity magnitude
    void pv_mag(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity in meridional plane
    void pv_merid(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// Computes the probability density distribution of the DM velocity along the azimuthal direction
    void pv_azim(int N, double R, double z, double *result, double tolerance=1e-3);
    
    /// GSL error handler
    static void GSL_error_func(const char * reason, const char * file, int line, int gsl_errno);

};
