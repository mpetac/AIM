#include "InversionInterp.hpp"

InversionInterp::InversionInterp(Model* model, double E, double Lz, double psiEnv, double h, int n) {
    std::complex<double> z2(1e4, 1e4);
    double theta[n], z2re[n], z2im[n];
    for (int i = 0; i < n; i++) {
        double t = M_PI * i / (n - 1.);
        std::complex<double> xi = 0.5 * psiEnv * (1. + std::cos(t) + 2.0i * h * sin(t));
        
        z2 = model->psi_inverse(xi, E, Lz, z2);
        theta[i] = t;
        z2re[i] = std::real(z2);
        z2im[i] = std::imag(z2);
    }
    
    reAcc = gsl_interp_accel_alloc();
    imAcc = gsl_interp_accel_alloc();
    z2Real = gsl_spline_alloc(gsl_interp_linear, n);
    z2Imag = gsl_spline_alloc(gsl_interp_linear, n);
    gsl_spline_init(z2Real, theta, z2re, n);
    gsl_spline_init(z2Imag, theta, z2im, n);
}


std::complex<double> InversionInterp::z2_eval(double t) {
    return gsl_spline_eval(InversionInterp::z2Real, t, InversionInterp::reAcc) + 1.0i * gsl_spline_eval(InversionInterp::z2Imag, t, InversionInterp::imAcc);
}

InversionInterp::~InversionInterp() {
    gsl_spline_free(z2Real);
    gsl_spline_free(z2Imag);
    gsl_interp_accel_free(reAcc);
    gsl_interp_accel_free(imAcc);
}
