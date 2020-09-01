#pragma once

#include <complex>
#include <gsl/gsl_spline.h>
#include "Model.hpp"

class InversionInterp {

private:
    gsl_spline *z2Real, *z2Imag;
    gsl_interp_accel *reAcc, *imAcc;
    
    Model *model;
    
public:
    InversionInterp(Model *model, double E, double Lz, double psiEnv, double h, int n);
    
    ~InversionInterp();
    
    std::complex<double> z2_eval(double t);
    
};
