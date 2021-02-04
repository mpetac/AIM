#pragma once

#include <time.h>
#include <future>
#include <vector>
#include <complex>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>

#include "InversionInterp.hpp"
#include "Model.hpp"
#include "halos/Halo.hpp"
#include "baryons/Baryons.hpp"

/**
 * Class used for computing the phase-space distribution function (PSDF).
 */ 

class Inversion {
    
private:
    /// Interpolation spline for the phase-space distribution function
    gsl_spline2d *F;
    /// Interpolation spline for the inverse of the angular momentum of circular orbits
    gsl_spline *LcInv;
    /// Accelerators for the evaluation of spline
    gsl_interp_accel *EAcc, *LzAcc, *LcInvAcc;
    /// Number of interpolation points
    size_t nInterpZ2 = 200., nInterpLc = 10000;
    /// Number of intervals used in the numerical integration
    size_t nIntervals = 1e5;
    
    /// Parameter controling verbose output
    bool verbose;
    /// Relative tolerance used in performing the contour integrals
    double tolerance_F;
    /// Width of the contours
    double h = 0.05;
    /// Normalization factor for the Lz-even part of PSDF
    double result_fact_even = 1. / (2. * pow(M_PI, 2) * std::sqrt(2.));
    /// Normalization factor for the Lz-odd part of PSDF
    double result_fact_odd = 1. / (8. * pow(M_PI, 2));
    
    /// Galactic model
    Model *model;
    /// Value of the total gravitational potential at the origin (i.e. its maximum value)
    double psi0;
    
    /// Computes the tabulated values of the PSDF
    void tabulate_F(int N_E, int N_Lz, double *Epts, double *Lzpts, double *Fpts);
    
    /// Computes the Lz-even part of the PSDF
    double F_even(double *params);
    
    /// Computes the Lz-odd part of the PSDF
    double F_odd(double *params);
    
    /// Computes PSDF using Eddington's inversion
    double F_eddington(double E);
    
    /// GSL error handler
    static void GSL_error_func(const char * reason, const char * file, int line, int gsl_errno);
    /// Silent GSL error handler
    static void GSL_error_func_silent(const char * reason, const char * file, int line, int gsl_errno);
    
public:
    /// Initializer which performs the interpolation of the PSDF
    Inversion(Model *model, int N_E, int N_Lz, double tolerance_F = 1e-3, bool verbose = 0);
    
    /// Destructor
    ~Inversion();
    
    /// Returns the value of the PSDF
    double eval_F(double E, double Lz);
    
    /// Returns the value of inverse circular velocity
    double eval_LcInv(double E);
};
