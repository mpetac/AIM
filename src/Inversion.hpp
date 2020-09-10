#pragma once

#include <thread>
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
    /// Interpolation spline for the maximum angular momentum at given relative energy
    gsl_spline *LcInv;
    /// Interpolation spline for the phase-space distribution function
    gsl_spline2d *F;
    /// Accelerators for the evaluation of splines
    gsl_interp_accel *EAcc, *LzAcc, *LcAcc;
    /// Number of intervals used in the numerical integration
    size_t nIntervals = 1e5;
    
    /// Parameter controling verbose output
    bool verbose;
    /// Number of interpolation points used in the inversion of total gravitational potential
    int nInterp = 200;
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
    
public:
    /// Initializer which performs the interpolation of the PSDF
    Inversion(Model *model, int N_E, int N_Lz, int N_Lc = 930, double tolerance_F = 1e-2, bool verbose = 0);
    
    /// Destructor
    ~Inversion();
    
    /// Returns the value of the PSDF
    double eval_F(double E, double Lz);
    
    /// Returns the inverse of maximum circular velocity
    double eval_LcI(double E);
    
    /// GSL error handler
    static void GSL_error_func(const char * reason, const char * file, int line, int gsl_errno);
    
    void test();
};
