#pragma once

#include <complex>
#include <gsl/gsl_spline.h>
#include "Model.hpp"

/**
 * Class used for interpolating the inverse of the total gravitational potential (i.e. z2 = psi^{-1}(xi, R2)
 */

class InversionInterp {

private:
    /// Splines for the real and imaginary values of z2
    gsl_spline *z2Real, *z2Imag;
    /// Accelerators for the evaluation of splines
    gsl_interp_accel *reAcc, *imAcc;
    
    /// Galactic model
    Model *model;
    
public:
    /// Initializer that interpolates the inverse of the total gravitational potential
    InversionInterp(Model *model, double E, double Lz, double psiEnv, double h, int n);
    
    /// Destructor
    ~InversionInterp();
    
    /// Provides the z2 value at a given point along the contour
    std::complex<double> z2_eval(double t);
    
};
