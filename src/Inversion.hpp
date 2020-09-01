#pragma once

#include <thread>
#include <future>
#include <vector>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include "InversionInterp.hpp"
#include "Model.hpp"
#include "halos/Halo.hpp"
#include "baryons/Baryons.hpp"

class Inversion {

private:
    gsl_spline *LcInv;
    gsl_spline2d *F;
    gsl_interp_accel *EAcc, *LzAcc, *LcAcc;
    size_t nIntervals = 1e5;
    
    bool verbose;
    int nInterp = 100;
    double tolerance_F;
    double h = 0.05;
    double result_fact_even = 1. / (2. * pow(M_PI, 2) * std::sqrt(2.));
    double result_fact_odd = 1. / (8. * pow(M_PI, 2));
    
    Model *model;
    double psi0;
    
    void tabulate_F(int N_E, int N_Lz, double *Epts, double *Lzpts, double *Fpts);
    double F_even(double *params);
    double F_odd(double *params);
    
public:
    Inversion(Model *model, int N_E, int N_Lz, int N_Lc = 930, double tolerance_F = 1e-3, bool verbose = 0);
    
    ~Inversion();
    
    static void GSL_error_func(const char * reason, const char * file, int line, int gsl_errno);
    
};
