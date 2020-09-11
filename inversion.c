#define _GNU_SOURCE 
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <complex.h>
#include <pthread.h> 
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>


//gcc inversion.c -o inversion.o -lgsl -lgslcblas -lm


/* --------------------------------
* Declerations of global variables
* -------------------------------- */

// Gravitational constant (set to 1, since I assume the input qunatities already carry the correct units)
const double G = 4.3e-6;

// Constants related to numerical integration & minimization
const double tolPsi = 1e-6, tolInverse = 1e-8, tolPSDF = 1e-2, tolEta = 1e-2;

// Constant that defines the thickness of the contour (its value is set later, depending on the particular relative energy and angular momentum)
double h = 101.769;

// Parameters related to motion of Sun and Earth
//double Usol = 0., Vsol = 0., Wsol = 0., LSR = 0., vE = 0., sE = 0, cE = 1;
double Usol = 11., Vsol = 12., Wsol = 7., LSR = 230., vE = 30., sE = 0, cE = 1;
double Rsol = 8.122;

// DM density profile parameters
double rhoS = 1., rS = 16., q = 1., alpha = 1., beta = 3., gama = 1., sin_e = 0;
double complex e = 0, casin_e = 0;

// Normalization of the spherioidal DM gravitational potential (avoid recomputing numerical integral for q != 1)
double psiDM0 = 0, psiDM0_num = 0;

// Constant that defines the fraction of Myamoto-Nagai disc potential (the rest of disc mass is modeled as and additional spherical Hernquist component)
double A_axi_thin = 1., a_axi_thin = 1., b_axi_thin = 0.1;
double A_axi_thick = 10., a_axi_thick = 10., b_axi_thick = 1.;
double A_thin = 1, B_thin = 0.1, A_thick = 1, B_thick = 0.1;
double k_A[5] = {-1.223844769510286715e-02, 1.475306058525469588e-01, -6.939601286064533969e-01, 2.682422756460722812e+00, -4.992464980932310049e-04};
double k_B[5] = {1.271287036483978738e+01, -3.330747680439659320e+01, 3.259178578793687819e+01, -1.456278256476855226e+01, 9.496121056454400744e+00};

// Hernquist gravitational potential parameters
double A_sph = 10., a_sph = 0.1;

// Constants related to the halo rotation
double vPhi = 1, rPhi2 = 1;

// Velocity power in astrophysical factors
double vPow = -1;

// Some additional global variables and functions...
int verbose = 0, tabulated = 0, core = 0;
clock_t lapse;
gsl_spline *z2_re, *z2_im;
gsl_interp_accel *z2_re_acc, *z2_im_acc;
gsl_spline2d *psdf;
gsl_spline *LcInv;
gsl_interp_accel *xAcc, *yAcc, *lAcc;
gsl_error_handler_t *errHand;


void (*eval_PSDF)(int, int, double*, double*, double*);
double (*Sv)(double, double, double, double);
double complex (*psi_exp_num_kernel)(double complex, double complex, double, double);

double complex (*rho)(double complex);
double complex (*drho_dm2)(double complex);
double complex (*drho_dz2)(double complex, double complex);
double complex (*d2rho_dz22)(double complex, double complex);

double complex (*psiDM)(double complex, double complex, double complex);
double complex (*dpsiDM_dR2)(double complex, double complex, double complex);
double complex (*dpsiDM_dz2)(double complex, double complex, double complex);
double complex (*d2psiDM_dz22)(double complex, double complex, double complex);
double complex (*d2psiDM_dz2dR2)(double complex, double complex, double complex);

struct stat st = {0};

// Sturcts used in passing the parameters in nested numerical integrals
struct psdf_params {
    double E;
    double L;
    double Rc;
    double psiEnv;
};

struct psi_q1_params {
    double complex r;
};

struct psi_params {
    double complex R2;
    double complex z2;
};

struct psi_exp_params {
    double complex R0;
    double complex z0;
    double RdInv;
    double zdInv;
    double R;
};

struct sigma_params {
    double R2;
};

struct params_1 {
    double P1;
};

struct params_9 {
    double P1;
    double P2;
    double P3;
    double P4;
    double P5;
    double P6;
    double P7;
    double P8;
    double P9;
};



/* --------------------------------------------------------
* GSL error handler (when enabled may throw many warnings) 
* -------------------------------------------------------- */

void errorHandlerFunc(const char * reason, const char * file, int line, int gsl_errno) {
    if (verbose) printf(" -> GSL error #%i in line %d of %s: %s\n", gsl_errno, line, file, gsl_strerror(gsl_errno));
}

double besselK0(double x) {
    if (x > 2e1) return 0;
    else return gsl_sf_bessel_K0(x);
}

double factorial(double n) {
    double result = 1;
    for (int i = 1; i < n + 1; i++) {
        result *= i;
    }
    return result;
}

double complex Ei(double complex x) {
    if (cimag(x) == 0) return gsl_sf_expint_Ei(x);
    else {
        double complex result = M_EULER + clog(x);
        double complex sum1 = 0;
        for (int i = 1; i < 100; i++) {
            double sum2 = 0;
            for (int j = 0; j < 0.5 * (i - 1); j++) {
                sum2 += 1. / (2 * j + 1);
            }
            sum1 += cpow(-1, i - 1) * cpow(x, i) * pow(2, 1 - i) / factorial(i) * sum2;
        }
        result += cexp(x / 2.) * sum1;
        return result;
    }
}


/* ----------------------------------------------------
* Definition of DM density profile and its derivatives 
* ---------------------------------------------------- */

double complex rhoNFW(double complex m2) {
    double complex x2 = m2 / (rS * rS);
    return rhoS * cpow(x2, -gama / 2.) * cpow(1. + cpow(x2, alpha / 2.), -(beta - gama) / alpha);
}

double complex drhoNFW_dm2(double complex m2) {
    double complex x = csqrt(m2) / rS;
    return -rhoNFW(m2) * (gama + beta * cpow(x, alpha)) / (2. * cpow(x, 2) * (1. + cpow(x, alpha))) * cpow(rS, -2); 
}

double complex drhoNFW_dz2(double complex R2, double complex z2) {
    double complex m2 = R2 + z2 * pow(q, -2);
    return drhoNFW_dm2(m2) * pow(q, -2);
}

double complex d2rhoNFW_dz22(double complex R2, double complex z2) {
    double complex m2 = R2 + z2 * pow(q, -2.);
    double complex m = csqrt(m2);
    double complex x = m / rS;
    return 0.25 * rhoNFW(m2) * cpow(m, -4) * cpow(1. + cpow(x, alpha), -2) * (gama * (2. + gama) + ((2. + alpha) * gama + beta * (2 - alpha + 2 * gama)) * cpow(x, alpha) + beta * (2 + beta) * cpow(x, 2 * alpha)) * pow(q, -4.);
}

double complex rhoBUR(double complex m2) {
    double complex x = csqrt(m2) / rS;
    return rhoS / ((1 + x) * (1 + cpow(x, 2)));
}

double complex drhoBUR_dm2(double complex m2) {
    double complex x = csqrt(m2) / rS;
    return - 0.5 * rhoS * (1 + 2 * x + 3 * cpow(x, 2)) / (x * cpow(1 + x + cpow(x, 2) + cpow(x, 3), 2)) * pow(rS, -2);
}

double complex drhoBUR_dz2(double complex R2, double complex z2) {
    double complex x = csqrt(R2 + z2 * pow(q, -2)) / rS;
    return - 0.5 * rhoS * (1 + 2 * x + 3 * cpow(x, 2)) / (x * cpow(1 + x + cpow(x, 2) + cpow(x, 3), 2)) * pow(rS * q, -2);
}

double complex d2rhoBUR_dz22(double complex R2, double complex z2) {
    double complex m2 = R2 + z2 * pow(q, -2.);
    double complex m = csqrt(m2);
    double complex x = m / rS;
    return 0.25 * rhoS * (1 + 3 * x * (1 + x * (2 + x * (6 + x * (7 + 5 * x))))) * cpow(x * (1 + x + cpow(x, 2) + cpow(x, 3)), -3) * pow(rS * q, -4);
}


/* -------------------------------------------------------------------------------------------
* Functions related to numerical integration of the DM gravitational potnetial in case q != 1
* ------------------------------------------------------------------------------------------- */

double psiDM0_integrand(double m2, void * params) {
    return creal(rho(m2));
}

double psiDM_integrand_q1_re(double t, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = (R2 + z2) * t;
    double complex V = (R2 + z2) * t * t;
    return creal(rho(U) * V * (1 - sqrt(t)) * pow(t, -2));
}

double psiDM_integrand_q1_im(double t, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = (R2 + z2) * t;
    double complex V = (R2 + z2) * t * t;
    return cimag(rho(U) * V * (1 - sqrt(t)) * pow(t, -2));
}

double psiDM_integrand_re(double t, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * t + z2 * cpow(q * q - 1. + 1. / t, -1);
    double complex V = R2 * t * t + z2 * cpow(q * q - 1. + 1. / t, -2);
    return creal(rho(U) * V * (casin_e - casin(e * sqrt(t))) * pow(t, -2));
}

double psiDM_integrand_im(double t, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * t + z2 * cpow(q * q - 1. + 1. / t, -1);
    double complex V = R2 * t * t + z2 * cpow(q * q - 1. + 1. / t, -2);
    return cimag(rho(U) * V * (casin_e - casin(e * sqrt(t))) * pow(t, -2));
}

double complex psiDM_num(double complex R2, double complex z2, double complex r) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    gsl_integration_workspace *workspace;
    gsl_function F;
    
    if (psiDM0_num == 0) {
        double result, abserr;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        F.function = &psiDM0_integrand;
        gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        psiDM0_num = 2. * M_PI * G * q * sin_e * result;
// 		printf("Psi_0: %g (%g +/- %g)\n", psiDM0_num, result, abserr);
    }
    
    if (R2 == 0 && z2 == 0) return psiDM0_num;
    else if (cabs(R2 + z2) > 1e16) return 0 + 0*I;
    else {
        double result_re, result_im, abserr_re, abserr_im;
        struct psi_params params = {R2, z2};
        F.params = &params;
        
        if (q == 1) F.function = &psiDM_integrand_q1_re;
        else F.function = &psiDM_integrand_re;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qag(&F, 0, 1, 0, tol, nIntervals, 6, workspace, &result_re, &abserr_re);
        gsl_integration_workspace_free(workspace);
        //printf("Psi_re: %g +/- %g\n", result_re, abserr_re);
        
        if (cimag(R2) != 0 || cimag(z2) != 0) {
            if (q == 1) F.function = &psiDM_integrand_q1_im;
            else F.function = &psiDM_integrand_im;
            workspace = gsl_integration_workspace_alloc(nIntervals);
            gsl_integration_qag(&F, 0, 1, 0, tol, nIntervals, 1, workspace, &result_im, &abserr_im);
            gsl_integration_workspace_free(workspace);
            //printf("Psi_im: %g +/- %g\n", result_im, abserr_im);
        } else result_im = 0;
        
        //printf("norm: %g +/- %g : re: %g +/- %g : im: %g +/- %g\n", result, abserr, result_re, abserr_re, result_im, abserr_im);
        
        return psiDM0_num - 2. * M_PI * G * q / e * (result_re + I * result_im);
    }
}

double dpsiDM_dR2_integrand_re(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return creal(rho(U) * cpow(1. + u, -2) * cpow(q * q + u, -0.5));
}

double dpsiDM_dR2_integrand_im(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return cimag(rho(U) * cpow(1. + u, -2) * cpow(q * q + u, -0.5));
}

double complex dpsiDM_dR2_num(double complex R2, double complex z2, double complex r) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result_re, result_im, abserr_re, abserr_im;
    struct psi_params params = {R2, z2};
    gsl_integration_workspace *workspace;
    gsl_function F;
    F.params = &params;
    
    F.function = &dpsiDM_dR2_integrand_re;
    workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace);
    
    if (cimag(R2) != 0 || cimag(z2) != 0) {
        F.function = &dpsiDM_dR2_integrand_im;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace);
    } else result_im = 0;
    
    return - M_PI * G * q * (result_re + I * result_im);
}

double dpsiDM_dz2_integrand_re(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return creal(rho(U) / (1. + u) * cpow(q * q + u, -3./2.));
}

double dpsiDM_dz2_integrand_im(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return cimag(rho(U) / (1. + u) * cpow(q * q + u, -3./2.));
}

double complex dpsiDM_dz2_num(double complex R2, double complex z2, double complex r) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result_re, result_im, abserr_re, abserr_im;
    struct psi_params params = {R2, z2};
    gsl_integration_workspace *workspace;
    gsl_function F;
    F.params = &params;
    
    F.function = &dpsiDM_dz2_integrand_re;
    workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace);
    
    if (cimag(R2) != 0 || cimag(z2) != 0) {
        F.function = &dpsiDM_dz2_integrand_im;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace);
    } else result_im = 0;
    
    return - M_PI * G * q * (result_re + I * result_im);
}

double d2psiDM_dz22_integrand_re(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return creal(drho_dm2(U) / (1. + u) * cpow(q * q + u, -5./2.));
}

double d2psiDM_dz22_integrand_im(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return cimag(drho_dm2(U) / (1. + u) * cpow(q * q + u, -5./2.));
}

double complex d2psiDM_dz22_num(double complex R2, double complex z2, double complex r) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result_re, result_im, abserr_re, abserr_im;
    struct psi_params params = {R2, z2};
    gsl_integration_workspace *workspace;
    gsl_function F;
    F.params = &params;
    
    F.function = &d2psiDM_dz22_integrand_re;
    workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace);
    
    if (cimag(R2) != 0 || cimag(z2) != 0) {
        F.function = &d2psiDM_dz22_integrand_im;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace);
    } else result_im = 0;
    
    return - M_PI * G * q * (result_re + I * result_im);
}

double d2psiDM_dz2dR2_integrand_re(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return creal(drho_dm2(U) * pow(1. + u, -2) * cpow(q * q + u, -3./2.));
}

double d2psiDM_dz2dR2_integrand_im(double u, void * params) {
    struct psi_params * p = (struct psi_params *) params;
    double complex R2 = (p->R2), z2 = (p->z2);
    double complex U = R2 * cpow(1. + u, -1) + z2 * cpow(q * q + u, -1);
    return creal(drho_dm2(U) * pow(1. + u, -2) * cpow(q * q + u, -3./2.));
}

double complex d2psiDM_dz2dR2_num(double complex R2, double complex z2, double complex r) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result_re, result_im, abserr_re, abserr_im;
    struct psi_params params = {R2, z2};
    gsl_integration_workspace *workspace;
    gsl_function F;
    F.params = &params;
    
    F.function = &d2psiDM_dz2dR2_integrand_re;
    workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace);
    
    if (cimag(R2) != 0 || cimag(z2) != 0) {
        F.function = &d2psiDM_dz2dR2_integrand_im;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace);
    } else result_im = 0;
    
    return - M_PI * G * q * (result_re + I * result_im);
}



/* ---------------------------------------------------------------------
* Functions related to numerical integration exponential disk potential
* --------------------------------------------------------------------- */

//psi
double complex psi_exp_kernel_0(double complex R0, double complex z0, double R, double z) {
    double complex arg = 2. * R / (csqrt(cpow(z - z0, 2) + cpow(R - R0, 2)) + csqrt(cpow(z - z0, 2) + cpow(R + R0, 2)));
    return casin(arg);
}

//dpsi_dR2
double complex psi_exp_kernel_1(double complex R0, double complex z0, double R, double z) {
    double complex Ap = csqrt(cpow(z - z0, 2) + cpow(R + R0, 2));
    double complex Am = csqrt(cpow(z - z0, 2) + cpow(R - R0, 2));
    return - R * ((R0 + R) / Ap + (R0 - R) / Am) / (R0 * csqrt(cpow(Ap + Am, 2) - 4. * pow(R, 2)) * (Ap + Am));
}

//dpsi_dz2
double complex psi_exp_kernel_2(double complex R0, double complex z0, double R, double z) {
    double complex Ap = csqrt(cpow(z - z0, 2) + cpow(R + R0, 2));
    double complex Am = csqrt(cpow(z - z0, 2) + cpow(R - R0, 2));
    return R * (z - z0) / (z0 * csqrt(cpow(Ap + Am, 2) - 4. * pow(R, 2)) * Ap * Am);
}

//d2psi_dz22
double complex psi_exp_kernel_3(double complex R0, double complex z0, double R, double z) {
    double complex Ap = csqrt(cpow(z - z0, 2) + cpow(R + R0, 2)), Am = csqrt(cpow(z - z0, 2) + cpow(R - R0, 2));
    double complex ApInv = 1. / Ap, AmInv = 1. / Am;
    double complex AsumInv = ApInv + AmInv;
    double complex AsumInv2 = cpow(AsumInv, 2);
    double complex C = (cpow(Ap + Am, 2) - 4 * pow(R, 2));
    double complex z2 = cpow(z0 - z, 2);
    return R * (-(Ap + Am) * C * (AsumInv - z2 * (1. / cpow(Ap, 3) + 1. / cpow(Am, 3))) + 4 * pow(R, 2) * z2 * AsumInv2 + 2 * z2 * AsumInv2 * C) / (2 * cpow(z0, 2) * cpow(C, 1.5) * cpow(Ap + Am, 3.5));
}

//d2psi_dz2dR2
double complex psi_exp_kernel_4(double complex R0, double complex z0, double R, double z) {
    double complex Ap = csqrt(cpow(z - z0, 2) + cpow(R + R0, 2)), Am = csqrt(cpow(z - z0, 2) + cpow(R - R0, 2));
    double complex ApInv = 1. / Ap, AmInv = 1. / Am;
    double complex AsumInv = ApInv + AmInv;
    double complex B = ((R0 + R) * ApInv + (R0 - R) * AmInv);
    double complex C = (cpow(Ap + Am, 2) - 4 * pow(R, 2));
    return R * (z0 - z) * (2 * pow(R, 2) * AsumInv * B + AsumInv * B * C + 0.5 * ((R0 + R) * cpow(Ap, 3) + (R0 - R) * cpow(Am, 3) * (Ap + Am) * C)) / (R0 * z0 * cpow(C, 1.5) * cpow(Ap + Am, 3.5));
}

double psi_exp_integrand2_re(double z, void * params) {
    struct psi_exp_params * p = (struct psi_exp_params *) params;
    double complex R0 = (p->R0), z0 = (p->z0);
    double RdInv = (p->RdInv), zdInv = (p->zdInv), R = (p->R);
    return creal(exp(-cabs(z) * zdInv) * psi_exp_num_kernel(R0, z0, R, z));
}

double psi_exp_num_re(void * params) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result, abserr;
    
    struct psi_exp_params * p = (struct psi_exp_params *) params;
    
    gsl_function F;
    F.params = params;
    F.function = &psi_exp_integrand2_re;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    if(psi_exp_num_kernel == psi_exp_kernel_2) gsl_integration_qags(&F, -1e2, 1e2, 0, tol, nIntervals, workspace, &result, &abserr);
    else gsl_integration_qagi(&F, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return result;
}

double psi_exp_integrand1_re(double R, void * params) {
    struct psi_exp_params * p = (struct psi_exp_params *) params;
    p->R = R;
    return R * besselK0(R * (p->RdInv)) * psi_exp_num_re(params);
}

double psi_exp_integrand2_im(double z, void * params) {
    struct psi_exp_params * p = (struct psi_exp_params *) params;
    double complex R0 = (p->R0), z0 = (p->z0);
    double RdInv = (p->RdInv), zdInv = (p->zdInv), R = (p->R);
    return cimag(exp(-cabs(z) * zdInv) * psi_exp_num_kernel(R0, z0, R, z));
}

double psi_exp_num_im(void * params) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result, abserr;
    
    gsl_function F;
    F.params = params;
    F.function = &psi_exp_integrand2_im;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qagi(&F, 0, tol, nIntervals, workspace, &result, &abserr);
    //gsl_integration_qags(&F, -1e8, 1e8, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return result;
}

double psi_exp_integrand1_im(double R, void * params) {
    struct psi_exp_params * p = (struct psi_exp_params *) params;
    p->R = R;
    return R * besselK0(R * (p->RdInv)) * psi_exp_num_im(params);
}

double complex psi_exp_num(double complex R2, double complex z2, double Rd, double zd, int derivative) {
    size_t nIntervals = 1e5;
    double tol = tolPsi;
    double result_re, result_im, abserr_re, abserr_im;
    struct psi_exp_params params = {csqrt(R2), csqrt(z2), 1. / Rd, 1. / zd, 0};
    
    switch(derivative) {
        case 0:
            psi_exp_num_kernel = &psi_exp_kernel_0;
            break;
        case 1:
            psi_exp_num_kernel = &psi_exp_kernel_1;
            break;
        case 2:
            psi_exp_num_kernel = &psi_exp_kernel_2;
            break;
        case 3:
            psi_exp_num_kernel = &psi_exp_kernel_3;
            break;
        case 4:
            psi_exp_num_kernel = &psi_exp_kernel_4;
            break;
    }
    
    gsl_integration_workspace *workspace;
    gsl_function F;
    F.params = &params;
    
    F.function = &psi_exp_integrand1_re;
    workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_re, &abserr_re);
    //gsl_integration_qags(&F, 1e-3, 10, 0, tol, nIntervals, workspace, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace);
    
    if (cimag(R2) != 0 || cimag(z2) != 0) {
        F.function = &psi_exp_integrand1_im;
        workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qagiu(&F, 0, 0, tol, nIntervals, workspace, &result_im, &abserr_im);
        //gsl_integration_qags(&F, 1e-3, 10, 0, tol, nIntervals, workspace, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace);
    } else result_im = 0;
    
    return 2. * G / (M_PI * pow(Rd, 3)) * (result_re + I * result_im);
}

/* -----------------------------------------------------------------
* Declerations of total gravitational potential and its derivatives
* ----------------------------------------------------------------- */

double complex psiNFW(double complex R2, double complex z2, double complex r) {
    if(r == 0) return 4. * M_PI * G * rhoS * pow(rS, 2);
    else return 4. * M_PI * G * rhoS * cpow(rS, 3) * clog(1 + r / rS) / r;
}

double complex dpsiNFW_dr2(double complex R2, double complex z2, double complex r) {
    return 2. * M_PI * G * rhoS * cpow(rS / r, 3) * (r / (r + rS) - clog(1 + r / rS));
}

double complex d2psiNFW_dr22(double complex R2, double complex z2, double complex r) {
    return M_PI * G * rhoS * pow(rS, 3) * (3. * cpow(r + rS, 2.) * clog(1 + r / rS) - 4. * cpow(r, 2.) - 3. * rS * r) / (cpow(r, 5.) * cpow(rS + r, 2));
}

double complex psiBUR(double complex R2, double complex z2, double complex r) {
    if(r == 0) return pow(M_PI, 2) * G * rhoS * pow(rS, 2);
    else {
        double complex x = r / rS;
        return G * M_PI * rhoS * pow(rS, 2) / x * (-M_PI + 2 * (1 + x) * (catan(1. / x) + clog(1 + x)) + (1 - x) * clog(1 + cpow(x, 2)));
    }
}

double complex dpsiBUR_dr2(double complex R2, double complex z2, double complex r) {
    double complex x = r / rS;
    return 0.5 * G * M_PI * rhoS * (M_PI - 2 * catan(1. / x) - clog(1 + cpow(x, 2)) - 2 * clog(1 + x)) * cpow(x, -3);
}

double complex d2psiBUR_dr22(double complex R2, double complex z2, double complex r) {
    double complex x = r / rS;
    return 0.25 * G * M_PI * rhoS * pow(rS, -2) * cpow(x, -5) / ((1 + x) * (1 + cpow(x, 2))) * (3 * (1 + x) * (cpow(x, 2) + 1) * (2 * catan(1. / x) + clog(1 + cpow(x, 2)) + 2 * clog(1 + x)) - 3 * M_PI * (cpow(x, 3) + cpow(x, 2) + x + 1) - 4 * cpow(x, 3));
}

double complex psiHer(double complex r, double A_Her, double a_Her) {
    return G * A_Her / (r + a_Her);
}

double complex dpsiHer_dr2(double complex r, double A_Her, double a_Her) {
    return - G * A_Her / (2. * r * cpow(r + a_Her, 2));
}

double complex d2psiHer_dr22(double complex r, double A_Her, double a_Her) {
    return G * A_Her * (3 * r + a_Her) / (4. * cpow(r, 3) * cpow(r + a_Her, 3));
}

double complex psiMN(double complex R2, double complex z2, double A_MN, double a_MN, double b_MN) {
    return G * A_MN / csqrt(R2 + cpow(a_MN + csqrt(z2 + b_MN * b_MN), 2));
}

double complex dpsiMN_dR2(double complex R2, double complex z2, double A_MN, double a_MN, double b_MN) {
    return - G * A_MN / (2. * cpow(R2 + cpow(a_MN + csqrt(z2 + b_MN * b_MN), 2), 3./2.));
}

double complex dpsiMN_dz2(double complex R2, double complex z2, double A_MN, double a_MN, double b_MN) {
    return - G * A_MN * (a_MN + csqrt(b_MN * b_MN + z2)) / (2. * csqrt(b_MN * b_MN + z2) * cpow(R2 + cpow(a_MN + csqrt(z2 + b_MN * b_MN), 2), 3./2.));
}

double complex d2psiMN_dz22(double complex R2, double complex z2, double A_MN, double a_MN, double b_MN) {
    return G * A_MN * (pow(a_MN , 3) + 5. * pow(a_MN , 2) * csqrt(b_MN * b_MN + z2) + 3. * cpow(b_MN * b_MN + z2, 3./2.) + a_MN * (7 * b_MN * b_MN + R2 + 7. * z2)) / (4. * cpow(b_MN * b_MN + z2, 3./2.) * cpow(R2 + cpow(a_MN + csqrt(z2 + b_MN * b_MN), 2), 5./2.));
}

double complex d2psiMN_dz2dR2(double complex R2, double complex z2, double A_MN, double a_MN, double b_MN) {
    return 3. * G * A_MN * (a_MN + csqrt(b_MN * b_MN + z2)) / (4. * csqrt(b_MN * b_MN + z2) * cpow(R2 + cpow(a_MN + csqrt(z2 + b_MN * b_MN), 2), 5./2.));
}

double complex psiExp(double complex x, double complex y) {	
    if (cabs(x) < 1e-1) x = 1e-1;
    if (cabs(y) < 1e-1) y = 1e-1;
    double complex val = 2. * G / M_PI * cexp(- x - y) * (M_EULER + Ei(x) - clog(x)) * (1 - M_EULER + Ei(y) - clog(y));
// 	if (isnan(creal(val))) printf("psiExp re nan at x = %g+%gI, y = %g+%gI -> %g %g %g \n", creal(x), cimag(x), creal(y), cimag(y), creal(cexp(- x - y)), creal(M_EULER + Ei(x) - clog(x)), creal(1 - M_EULER + Ei(y) - clog(y)));
// 	if (isnan(cimag(val))) printf("psiExp im nan at x = %g+%gI, y = %g+%gI -> %g %g %g \n", creal(x), cimag(x), creal(y), cimag(y), cimag(cexp(- x - y)), cimag(M_EULER + Ei(x) - clog(x)), cimag(1 - M_EULER + Ei(y) - clog(y)));
// 	if (isnan(creal(val))) printf("psiExp re nan at x = %g+%gI, y = %g+%gI -> %g %g %g \n", creal(x), cimag(x), creal(y), cimag(y), creal(cexp(-y)), creal(Ei(y)), creal(clog(y)));
// 	if (isnan(cimag(val))) printf("psiExp im nan at x = %g+%gI, y = %g+%gI -> %g %g %g \n", creal(x), cimag(x), creal(y), cimag(y), cimag(cexp(-y)), cimag(Ei(y)), cimag(clog(y)));
    return val;
}

double complex dpsiExp_dR2(double complex x, double complex y) {
    if (cabs(x) < 1e-8) x = 1e-8;
    if (cabs(y) < 1e-8) y = 1e-8;
    double complex val = - G / (M_PI * cpow(x, 2)) * cexp(- x - y) * (1 - cexp(x) + x * (M_EULER + Ei(x) - clog(x))) * (1 - M_EULER + Ei(y) - clog(y));
// 	if (isnan(creal(val))) printf("dpsiExp_dR2 nan at x, y: %g %g\n", creal(x), creal(y));
    return val;
}

double complex dpsiExp_dz2(double complex x, double complex y) {
    if (cabs(x) < 1e-8) x = 1e-8;
    if (cabs(y) < 1e-8) y = 1e-8;
    double complex val = - G / (M_PI * cpow(y, 2)) * cexp(- x - y) * (M_EULER + Ei(x) - clog(x)) * (1 - cexp(y) + y * (1 - M_EULER + Ei(y) - clog(y)));
// 	if (isnan(creal(val))) printf("dpsiExp_dz2 nan at x, y: %g %g\n", creal(x), creal(y));
    return val;
}

double complex d2psiExp_dz22(double complex x, double complex y) {
    if (cabs(x) < 1e-8) x = 1e-8;
    if (cabs(y) < 1e-8) y = 1e-8;
    double complex val = G / (2. * M_PI * cpow(y, 4)) * cexp(- x - y) * (M_EULER + Ei(x) - clog(x)) * ((1. + y) * (2. + y * (1. - M_EULER + Ei(y) - clog(y))) - (2. + y) * cexp(y));
// 	if (isnan(creal(val))) printf("d2psiExp_dz22 nan at x, y: %g %g\n", creal(x), creal(y));
    return val;
}

double complex d2psiExp_dz2dR2(double complex x, double complex y) {
    if (cabs(x) < 1e-8) x = 1e-8;
    if (cabs(y) < 1e-8) y = 1e-8;
    double complex val = G / (2. * M_PI * pow(x * y, 2)) * cexp(- x - y) * (1 - cexp(x) + x * (M_EULER + Ei(x) - clog(x))) * (1 - cexp(y) + y * (1 - M_EULER + Ei(y) - clog(y)));
// 	if (isnan(creal(val))) printf("d2psiExp_dz2dR2 nan at x, y: %g %g\n", creal(x), creal(y));
    return val;
}

double complex psiDisk(double complex R2, double complex z2) {
    return psiMN(R2, z2, A_axi_thin, a_axi_thin, b_axi_thin) + psiMN(R2, z2, A_axi_thick, a_axi_thick, b_axi_thick);
    
// 	double complex R = csqrt(R2), z = csqrt(z2);
// 	return A_axi_thin * A_thin * psiExp(R / a_axi_thin, z / b_axi_thin / B_thin) / a_axi_thin + A_axi_thick * A_thick * psiExp(R / a_axi_thick, z / b_axi_thick / B_thick) / a_axi_thick;
}

double complex dpsiDisk_dR2(double complex R2, double complex z2) {
    return dpsiMN_dR2(R2, z2, A_axi_thin, a_axi_thin, b_axi_thin) + dpsiMN_dR2(R2, z2, A_axi_thick, a_axi_thick, b_axi_thick);
    
// 	double complex R = csqrt(R2), z = csqrt(z2);
// 	return A_axi_thin * A_thin * dpsiExp_dR2(R / a_axi_thin, z / b_axi_thin / B_thin) * pow(a_axi_thin, -3) + A_axi_thick * A_thick * dpsiExp_dR2(R / a_axi_thick, z / b_axi_thick / B_thick) * pow(a_axi_thick, -3);
}

double complex dpsiDisk_dz2(double complex R2, double complex z2) {
    return dpsiMN_dz2(R2, z2, A_axi_thin, a_axi_thin, b_axi_thin) + dpsiMN_dz2(R2, z2, A_axi_thick, a_axi_thick, b_axi_thick);
    
// 	double complex R = csqrt(R2), z = csqrt(z2);
// 	return A_axi_thin * A_thin * dpsiExp_dz2(R / a_axi_thin, z / b_axi_thin / B_thin) / a_axi_thin * pow(b_axi_thin * B_thin, -2) + A_axi_thick * A_thick * dpsiExp_dz2(R / a_axi_thick, z / b_axi_thick / B_thick) / a_axi_thick * pow(b_axi_thick * B_thick, -2);
}

double complex d2psiDisk_dz22(double complex R2, double complex z2) {
    return d2psiMN_dz22(R2, z2, A_axi_thin, a_axi_thin, b_axi_thin) + d2psiMN_dz22(R2, z2, A_axi_thick, a_axi_thick, b_axi_thick);
    
// 	double complex R = csqrt(R2), z = csqrt(z2);
// 	return A_axi_thin * A_thin * d2psiExp_dz22(R / a_axi_thin, z / b_axi_thin / B_thin) / a_axi_thin * pow(b_axi_thin * B_thin, -4) + A_axi_thick * A_thick * d2psiExp_dz22(R / a_axi_thick, z / b_axi_thick / B_thick) / a_axi_thick * pow(b_axi_thick * B_thick, -4);
}

double complex d2psiDisk_dz2dR2(double complex R2, double complex z2) {
    return d2psiMN_dz2dR2(R2, z2, A_axi_thin, a_axi_thin, b_axi_thin) + d2psiMN_dz2dR2(R2, z2, A_axi_thick, a_axi_thick, b_axi_thick);
    
// 	double complex R = csqrt(R2), z = csqrt(z2);
// 	return A_axi_thin * A_thin * d2psiExp_dz2dR2(R / a_axi_thin, z / b_axi_thin / B_thin) / pow(a_axi_thin, -3) * pow(b_axi_thin * B_thin, -2) + A_axi_thick * A_thick * d2psiExp_dz2dR2(R / a_axi_thick, z / b_axi_thick / B_thick) * pow(a_axi_thick, -3) * pow(b_axi_thick * B_thick, -2);
}

double complex psi(double complex R2, double complex z2) {
    double complex r = csqrt(R2 + z2);
    return psiDM(R2, z2, r) + psiDisk(R2, z2) + psiHer(r, A_sph, a_sph);
}

double complex dpsi_dR2(double complex R2, double complex z2) {
    double complex r = csqrt(R2 + z2);
    return dpsiDM_dR2(R2, z2, r) + dpsiDisk_dR2(R2, z2) + dpsiHer_dr2(r, A_sph, a_sph);
}

double complex dpsi_dz2(double complex R2, double complex z2) {
    double complex r = csqrt(R2 + z2);
    return dpsiDM_dz2(R2, z2, r) + dpsiDisk_dz2(R2, z2) + dpsiHer_dr2(r, A_sph, a_sph);
}

double complex d2psi_dz22(double complex R2, double complex z2) {
    double complex r = csqrt(R2 + z2);
    return d2psiDM_dz22(R2, z2, r) + d2psiDisk_dz22(R2, z2) + d2psiHer_dr22(r, A_sph, a_sph);
}

double complex d2psi_dz2dR2(double complex R2, double complex z2) {
    double complex r = csqrt(R2 + z2);
    return d2psiDM_dz2dR2(R2, z2, r) + d2psiDisk_dz2dR2(R2, z2) + d2psiHer_dr22(r, A_sph, a_sph);
}



/* -------------------------------------------------------------------
* Functions used for computing the inverse of gravitational potential
* ------------------------------------------------------------------- */

double psi_inverse_func (const gsl_vector *v, void *params) {
    double complex *p = (double complex *)params;
    double complex z2 = gsl_vector_get(v, 0) + I * gsl_vector_get(v, 1);
    double complex psiVal = psi(p[1], z2);
    //printf("Psi(%g+%gI, %g+%gI): %g+%gI -> diff: %.32g\n", creal(p[1]), cimag(p[1]), creal(z2), cimag(z2), creal(psiVal), cimag(psiVal), cabs((p[0] - psiVal) / p[0]));
    return cabs((p[0] - psiVal) / p[0]);
}

double complex psi_inverse(double complex xi, double E, double L, double zr, double zi) {
    int iter = 0, lim = 200;
    double diff = 1, tol = tolInverse;
    //double complex par[3] = {xi, E, L};
    double complex par[2] = {xi, pow(L, 2) / (2. * (xi - E))};
    //printf("R2: %g+%gI\n", creal(par[1]), cimag(par[1]));
    
    // For large Re(z) and Im(z) use more iterations
    if(zr == 1e4 && zi == 1e4) lim = 1000;
    
    gsl_multimin_function my_func;
    my_func.n = 2;
    my_func.f = psi_inverse_func;
    my_func.params = par;

    gsl_vector *x, *step;
    x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, zr);
    gsl_vector_set(x, 1, zi);
    step = gsl_vector_alloc(2);
    gsl_vector_set(step, 0, fabs(zr + 1e-1) / 10.);
    gsl_vector_set(step, 1, fabs(zi + 1e-1) / 10.);
    //gsl_vector_set(step, 0, 0.1);
    //gsl_vector_set(step, 1, 0.1);

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 2);
    gsl_multimin_fminimizer_set(s, &my_func, x, step);

    while (diff > tol && iter < lim) {
        iter++;
        gsl_multimin_fminimizer_iterate(s);
        diff = s->fval;
    }
    
    double complex result = gsl_vector_get(s->x, 0) + gsl_vector_get (s->x, 1) * I;
    gsl_vector_free(x);
    gsl_vector_free(step);
    gsl_multimin_fminimizer_free(s);
    
    if (diff > tol) {
        //if(verbose) printf("Inversion failed! (xi: %.15g+%.15gI, E: %.15g, L: %.15g, iter: %d, delta: %g, z0: %g+%gI -> %g+%gI)\n", creal(xi), cimag(xi), E, L, iter, s->fval, zr, zi, creal(result), cimag(result));
        //gsl_vector *y;
        //y = gsl_vector_alloc(2);
        //gsl_vector_set(y, 0, creal(result));
        //gsl_vector_set(y, 1, cimag(result));
        //result = psi_inverse(xi, E, L, creal(result), cimag(result));
        //gsl_vector_free(y);
    } //else if(verbose) printf("Success in %d iterations!  (xi: %g+%gI, E: %g, L: %g, delta: %g)\n", iter, creal(xi), cimag(xi), E, L, s->fval);
    //printf("Result: %g+%gI (%d iterations)\n", creal(result), cimag(result), iter);
    return result;
}



/* ----------------------------------------------------------------------------
* Function used to compute the radius of circular orbit with relative energy E
* 
* Note: R_circ(E) is monotonic function
* ---------------------------------------------------------------------------- */

double Rcirc(double E) {
    int iter = 0, lim=1000;
    double Rc = 0., tol = tolInverse, dR = rS;
    while (iter < lim) {
        iter++;
        double R = Rc + dR;
        double Er = creal(psi(R * R, 0) + R * R * dpsi_dR2(R * R, 0));
        
        if (E > Er) dR /= 2.;
        else {
            Rc = R;
            dR *= 2;
        }
        
        //printf("E: %g, Er: %g -> delta: %g, Rc: %g, R: %g, dR: %g\n", E, Er, fabs(Er - E) / (tol * E), Rc, R, dR);
        if (fabs(Er - E) < tol) {
            Rc = R;
            break;
        }
        
    }
    if (verbose && iter >= lim) printf("Rc inversion did not converge! (E: %g)\n", E);
    //else if (verbose) printf("Rc inversion successful! (psiDM0: %g, E: %g, Rc: %g)\n", cimag(psiDM0), E, Rc);
    //else printf("N iter: %d\n", iter);
    return Rc;
}



/* ----------------------------------------------------------------------
* Second derivative of the DM density w.r.t. the gravitational potential
* ---------------------------------------------------------------------- */

double complex d2rho_dpsi2(double complex xi, double t, double E, double L) {
    double complex result = 0;
    if(cabs(xi) / creal(psiDM0) > 1e-5) {
        double complex R2 = pow(L, 2) / (2. * (xi - E));
        double complex z2;
        if (tabulated) z2 = psi_inverse(xi, E, L, gsl_spline_eval(z2_re, t, z2_re_acc), gsl_spline_eval(z2_im, t, z2_im_acc));
        else z2 = psi_inverse(xi, E, L, 1e3, 1e3);
        //else z2 = psi_inverse(xi, E, L, -8.5755e+37, -7.15295e+26);
        
        //printf("%g: %g+%gI -> %g+%gI\n", t, gsl_spline_eval(z2_re, t, z2_re_acc),  gsl_spline_eval(z2_im, t, z2_im_acc), creal(z2), cimag(z2));
        //printf("R2: %.8g+%.8gI, z2: %.8g+%.8gI\n", creal(R2), cimag(R2), creal(z2), cimag(z2));
        
        double complex dpsi = dpsi_dz2(R2, z2);
        result = d2rho_dz22(R2, z2) * cpow(dpsi, -2) - drho_dz2(R2, z2) * d2psi_dz22(R2, z2) * cpow(dpsi, -3);
    }
    return result;
}



/* ----------------------------------------------------------------------------------
* Tabulation of z^2 along the contour for given relative energy and angular momentum
* 
* Note: weather it uses ellipse or box conture is must be CONTROLED MANUALY by 
*       (un)commenting the corresponding lines!
* Note: for box contour it does not lead to significant speed-up
* ---------------------------------------------------------------------------------- */

void tabulate_z2(double E, double L, double psiEnv, int n) {
    clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    double complex z2 = 1e4 + 1e4 * I;
    double theta[n];
    double z2re[n];
    double z2im[n];
    for (int i = 0; i < n; i++) {
        double t = M_PI * i / (n - 1.);
        double complex xi = 0.5 * psiEnv * (1. + cos(t)) + I * h * sin(t);
        //double complex xi = psiEnv + I * h * i / (1. * n);
        //double complex xi = psiEnv * (n - i) / (1. * n) + I * h;
        //double complex xi = 0 + I * h * (n - i) / (1. * n);
        
        //printf("xi: %g + %g I\n", creal(xi), cimag(xi));
        
        z2 = psi_inverse(xi, E, L, creal(z2), cimag(z2));
        theta[i] = t;
        z2re[i] = creal(z2);
        z2im[i] = cimag(z2);
        
        //double complex R2 = L * L / (2. * (xi - E));
        //double complex w1 = z2 + b_axi*b_axi;
        //double complex w2 = R2 + cpow(a_axi + csqrt(z2 + b_axi*b_axi), 2);
        //printf("|z2|: %g, fi: %g\t|w1|: %g, fi1: %g\t|w2|: %g, fi2: %g\n", cabs(z2), carg(z2), cabs(w1), carg(w1), cabs(w2), carg(w2));
        
        //printf("|z2|: %g, fi: %g\n", cabs(z2), carg(z2));
        
    }
    
    //z2_re = gsl_spline_alloc(gsl_interp_cspline, n);
    //z2_im = gsl_spline_alloc(gsl_interp_cspline, n);
    z2_re = gsl_spline_alloc(gsl_interp_linear, n);
    z2_im = gsl_spline_alloc(gsl_interp_linear, n);
    gsl_spline_init(z2_re, theta, z2re, n);
    gsl_spline_init(z2_im, theta, z2im, n);
    z2_re_acc = gsl_interp_accel_alloc();
    z2_im_acc = gsl_interp_accel_alloc();
    
    /*
    for (int i = 0; i < n; i++) {
        double t = M_PI * i / (n - 1);
        //double t = M_PI * log(1. + i) / log(n);
        double complex xi = 0.5 * psiEnv * (1. + cos(t)) + I * h * sin(t);
        printf("z2(%g | %g+%gI): %g / %g, %g / %g\n", t, creal(xi), cimag(xi), z2re[i], gsl_spline_eval(z2_re, t, z2_re_acc), z2im[i], gsl_spline_eval(z2_im, t, z2_im_acc));
    }
    */
    tabulated = 1;
    if (verbose) printf("Tabulation done in %g s\n", (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
}



/* ------------------------------------------------------------------------
* Functions used to compute the L_z-even PSDF along the elliptical contour
* ------------------------------------------------------------------------ */

double PSDF_integrand(double t, void * params) {
    struct psdf_params * p = (struct psdf_params *) params;
    double E = (p->E), L = (p->L), Rc = (p->Rc), psiEnv = (p->psiEnv);
    double complex xi = 0.5 * psiEnv * (1. + cos(t)) + I * h * sin(t);
    double result = (0.5 * psiEnv * sin(t) * I + h * cos(t)) * cpow(xi - E, -0.5) * d2rho_dpsi2(xi, t, E, L);
    //printf("%g+%gI, %g+%gI, %g+%gI\n", creal((0.5 * psiEnv * sin(t) * I + h * cos(t))), cimag((0.5 * psiEnv * sin(t) * I + h * cos(t))), creal(cpow(xi - E, -0.5)), cimag(cpow(xi - E, -0.5)), creal(d2rho_dpsi2(xi, t, E, L)), cimag(d2rho_dpsi2(xi, t, E, L)));
    return creal(result);
}



double PSDF_even(double E, double L) {
    size_t nIntervals = 1e5;
    double tol = tolPSDF;
    double result, abserr;
    
    if (E == 0) return 0;
    
    E = E * psi(0, 0);
    double Rc = Rcirc(E);
    double psiEnv = psi(Rc * Rc, 0);
    L = L * pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0));
    h = 0.05 * psiEnv;
    
    //if (verbose) printf("Evaluating PSDF(%g, %g; %g)\n", E, L, psiDM0);
    
    tabulate_z2(E, L, psiEnv, 200);
    struct psdf_params params = {E, L, Rc, psiEnv};
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_function F;
    F.params = &params;
    F.function = &PSDF_integrand;
    gsl_integration_qags(&F, 0, M_PI, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    gsl_spline_free(z2_re);
    gsl_spline_free(z2_im);
    gsl_interp_accel_free(z2_re_acc);
    gsl_interp_accel_free(z2_im_acc);
    
    if (verbose && abserr / result > tol) printf("Warning: error larger then tolerance for PSDF(%g, %g)\n", E, L);
    return result / (2. * pow(M_PI, 2) * sqrt(2));
}



/* -----------------------------------------------------------------
* Functions used to compute the L_z-even PSDF along the box contour
* ----------------------------------------------------------------- */

double PSDF_integrand1(double t, void * params) {
    struct psdf_params * p = (struct psdf_params *) params;
    double E = (p->E), L = (p->L), Rc = (p->Rc), psiEnv = (p->psiEnv);
    double complex xi = psiEnv + I * t;
    double result = cpow(xi - E, -0.5) * d2rho_dpsi2(xi, t, E, L);
    return creal(result);
}

double PSDF_integrand2(double t, void * params) {
    struct psdf_params * p = (struct psdf_params *) params;
    double E = (p->E), L = (p->L), Rc = (p->Rc), psiEnv = (p->psiEnv);
    double complex xi = t + I * h;
    double result = I * cpow(xi - E, -0.5) * d2rho_dpsi2(xi, t, E, L);
    return creal(result);
}

double PSDF_integrand3(double t, void * params) {
    struct psdf_params * p = (struct psdf_params *) params;
    double E = (p->E), L = (p->L), Rc = (p->Rc), psiEnv = (p->psiEnv);
    double complex xi = I * t;
    double result = - cpow(xi - E, -0.5) * d2rho_dpsi2(xi, t, E, L);
    return creal(result);
}

double PSDF_box(double E, double L) {
    size_t nIntervals = 1e5;
    double tol = tolPSDF;
    double result1, abserr1, result2, abserr2, result3, abserr3;
    
    E = E * psi(0, 0);
    double Rc = Rcirc(E);
    double psiEnv = psi(Rc * Rc, 0);
    L = L * pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0));
    h = 0.05 * psiEnv;
    
    if (verbose) printf("Evaluating PSDF_box(%g, %g, %g %g)\n", E, L, psiEnv, Rc);
    
    tabulated = 0;
    //tabulate_z2(E, L, psiEnv, 200);
    
    
    struct psdf_params params = {E, L, Rc, psiEnv};
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_function F;
    F.params = &params;
    F.function = &PSDF_integrand1;
    gsl_integration_qags(&F, 0, h, 0, tol, nIntervals, workspace, &result1, &abserr1);
    gsl_integration_workspace_free(workspace);
    
    workspace = gsl_integration_workspace_alloc(nIntervals);
    F.function = &PSDF_integrand2;
    gsl_integration_qags(&F, 0, psiEnv, 0, tol, nIntervals, workspace, &result2, &abserr2);
    gsl_integration_workspace_free(workspace);
    
    
    result3 = 0;
    abserr3 = 0;
    //workspace = gsl_integration_workspace_alloc(nIntervals);
    //F.function = &PSDF_integrand3;
    //gsl_integration_qags(&F, 0, h, 0, tol, nIntervals, workspace, &result3, &abserr3);
    //gsl_integration_workspace_free(workspace);
    
    //printf("Integral 1: %g (1 +/- %g) done in %gs\n", result1, abserr1 / result1, (lapse2 - lapse1) / 1000.);
    //printf("Integral 2: %g (1 +/- %g) done in %gs\n", result2, abserr2 / result2, (lapse3 - lapse2) / 1000.);
    //printf("Integral 3: %g (1 +/- %g) done in %gs\n", result3, abserr3 / result3, (clock() / (CLOCKS_PER_SEC / 1000) - lapse3) / 1000.);
    
    return (result1 + result2 + result3) / (2. * pow(M_PI, 2) * sqrt(2));
}



/* -----------------------------------------------------------------
* Halo rotational velocity profile (that we used in our last paper)
* ----------------------------------------------------------------- */

double complex vPhiFunc(double complex xi, double E, double L) {
    return vPhi * rPhi2 * L * csqrt(2 * (xi - E)) / (2 * (xi - E) * rPhi2 + pow(L, 2));
//     return vPhi * clog(1. + cpow(pow(L, 2) / (2. * (xi - E) * rPhi2)), 2);
} 



/* -----------------------------------------------------------------------
* Functions used to compute the L_z-odd PSDF along the elliptical contour
* ----------------------------------------------------------------------- */

double PSDF_integrand_odd(double t, void * params) {
    struct psdf_params * p = (struct psdf_params *) params;
    double E = (p->E), L = (p->L), Rc = (p->Rc), psiEnv = (p->psiEnv);
    double complex xi = 0.5 * psiEnv * (1. + cos(t)) + I * h * sin(t);
    double result = (0.5 * psiEnv * sin(t) * I + h * cos(t)) / (xi - E) * d2rho_dpsi2(xi, t, E, L) * vPhiFunc(xi, E, L);
    return creal(result);
}

double PSDF_odd(double E, double L) {
    size_t nIntervals = 1e5;
    double tol = tolPSDF;
    double result, abserr;
    
    if (E == 0) return 0;
    
    E = E * psi(0, 0);
    double Rc = Rcirc(E);
    double psiEnv = psi(Rc * Rc, 0);
    L = L * pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0));
    h = 0.05 * psiEnv;
    if (verbose) printf("Evaluating PSDF(%g, %g, %g %g)\n", E, L, Rc, psiEnv);
    tabulate_z2(E, L, psiEnv, 200);
    
    struct psdf_params params = {E, L, Rc, psiEnv};
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_function F;
    F.params = &params;
    F.function = &PSDF_integrand_odd;
    gsl_integration_qags(&F, 0, M_PI, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    gsl_spline_free(z2_re);
    gsl_spline_free(z2_im);
    gsl_interp_accel_free(z2_re_acc);
    gsl_interp_accel_free(z2_im_acc);
    
    if (verbose && abserr / result > tol) printf("Warning: error larger then tolerance for PSDF(%g, %g)\n", E, L);
    return result / (8. * pow(M_PI, 2));
}

/*
void * sym_PSDF_thread(void *vargp) {
    double *args = (double *) vargp; 
    printf("E = %g, L = %g\n", args[0], args[1]);
    double *r = malloc(sizeof(double)); 
    double res = PSDF_even(args[0], args[1]);
    if (verbose) printf("E = %g, L = %g -> %g\n", args[0], args[1], res);
    *r = res;
    return (void *) r;
} 

void sym_PSDF(int n, int m, double *x, double *y, double *z) {
    int num = n * m;
    double *vals = malloc(num * sizeof(double));
    pthread_t thread[num];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double *args = malloc(2 * sizeof(double));
            args[0] = x[i];
            args[1] = y[m - 1 + j];
            pthread_create(&thread[i * m + j], NULL, sym_PSDF_thread, (void*)args);
            free(args);
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            void *result;
            printf("joining thread: %d\n", pthread_join(thread[i * m + j], (void **) &result));
            double val_e = *((double *) result);
            printf("val: %g\n", val_e);
            if (val_e < 0) val_e = 0;
            z[(m - 1 + j) * n + i] = log(1 + val_e);
            z[(m - 1 - j) * n + i] = log(1 + val_e);
            free(result);
        }
    }
    free(vals);
}
*/
/*
void * asym_PSDF_thread(void *vargp) {
    double *args = (double *) vargp; 
    double *r = malloc(sizeof(double)); 
    *r = PSDF_even(args[0], args[1]);
    if (verbose) printf("E = %g, L = %g -> %g\n", args[0], args[1], *r);
    return (void *) r;
} 

void asym_PSDF(int n, int m, double *x, double *y, double *z) {
    int num = n * m;
    double *vals = malloc(num * sizeof(double));
    pthread_t thread[num];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double *args = malloc(2 * sizeof(double));
            args[0] = x[i];
            args[1] = y[m - 1 + j];
            pthread_create(&thread[i * m + j], NULL, sym_PSDF_thread, (void*)args);
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            void *result;
            pthread_join(thread[i * m + j], (void **) &result);
            double val_e = * (double *) result;
            if (val_e < 0) val_e = 0;
            z[(m - 1 + j) * n + i] = log(1 + val_e);
            z[(m - 1 - j) * n + i] = log(1 + val_e);
        }
    }
}
*/


void sym_PSDF(int n, int m, double *x, double *y, double *z) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double val_e = PSDF_even(x[i], y[m - 1 + j]);
            if (val_e < 0) {
// 				printf("Negative PSDF: %g, %g, %g\n", val_e, x[i], y[m - 1 + j]);
                val_e = 0;
            }
                
            double val_o = 0;// vPhi * val_e * y[m - 1 + j];
            z[(m - 1 + j) * n + i] = log(1 + val_e + val_o);
            z[(m - 1 - j) * n + i] = log(1 + val_e - val_o);
            
            printf("PSDF computed: %g, %g -> %g, %g\n", x[i], y[m - 1 + j], val_e, val_o);
        }
    }
}
/**/


void asym_PSDF(int n, int m, double *x, double *y, double *z) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double val_e = PSDF_even(x[i], y[m - 1 + j]);
            double val_o = PSDF_odd(x[i], y[m - 1 + j]);
            if (val_e < 0) val_e = 0;
            if (val_o < 0) val_o = 0;
            if (val_o > val_e) val_o = val_e;
            z[(m - 1 + j) * n + i] = log(1 + val_e + val_o);
            z[(m - 1 - j) * n + i] = log(1 + val_e - val_o);
            
            printf("PSDF computed: %g, %g -> %g, %g\n", x[i], y[m - 1 + j], val_e, val_o);
        }
    }
}




/* --------------------------------------
* Functions used to interpolate the PSDF
* -------------------------------------- */

void interpolate_PSDF(int n, int m) {
    double x[n];
    double y[2 * m - 1];
    double z[n * (2 * m - 1)];
    
    psiDM0 = psi(0, 0);
    
    for (int i = 0; i < n; i++) {
        x[i] = pow(2 - 1e-4, 1. * i / (n - 1)) - 1;
    }
    
    for (int j = 0; j < 2 * m - 1; j++) {
        y[j] = 2. * j / (2. * m - 2.) - 1.;
    }
    
    eval_PSDF(n, m, x, y, z);
    
    psdf = gsl_spline2d_alloc(gsl_interp2d_bilinear, n, 2 * m - 1);
    gsl_spline2d_init(psdf, x, y, z, n, 2 * m - 1);
    xAcc = gsl_interp_accel_alloc();
    yAcc = gsl_interp_accel_alloc();
    
    double lx[5 * n];
    double l[5 * n];
    
    double r = 1.02;
    for (int i = 0; i < 5 * n; i++) {
        lx[5 * n - 1 - i] = 1 - pow(r, i) / pow(r, 5 * n - 1);
    }
    
    for (int i = 0; i < 5 * n; i++) {
        double Rc = Rcirc(lx[i] * psiDM0);
        l[i] = 1. / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
// 		printf("Lc(%g): %g\n", lx[i], 1. / l[i]);
    }
    LcInv = gsl_spline_alloc(gsl_interp_cspline, 5 * n);
    gsl_spline_init(LcInv, lx, l, 5 * n);
    lAcc = gsl_interp_accel_alloc();
    
    tabulated = 1;
    if (verbose) printf("Interpolation done!\n");
}




/* ----------------------------------------------------------------
* Functions used to compute the PSDF reconstructed density profile
* ----------------------------------------------------------------- */

double rho_num_c_integrand(double c, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = p->P1;
    double psiRz = p->P2;
    double v = p->P3;
    
    double E = (psiRz - 0.5 * pow(v, 2)) / psi(0, 0);
    
    //double L = R * v * c * gsl_spline_eval(LcInv, E, lAcc);
    double Rc = Rcirc(E * psi(0, 0));
    double L = R * v * c / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
    
    return exp(gsl_spline2d_eval(psdf, E, L, xAcc, yAcc)) - 1;
}

double rho_num_c_integral(void * params) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 * p = (struct params_9 *) params;
    
    gsl_function F;
    F.function = &rho_num_c_integrand;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, 1, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Phi int result: %g\n", result);
    return result;
}

double rho_num_v_integrand(double v, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    p->P3 = v;
    return pow(v, 2) * rho_num_c_integral(p);
}

double rho_num(double R2, double z2) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    double psiRz = creal(psi(R2, z2));
    
    struct params_9 p = {sqrt(R2), psiRz};
    
    gsl_function F;
    F.function = &rho_num_v_integrand;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, sqrt(2 * psiRz), 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return 4 * M_PI * result;
}

double PSDF_binned_vm_integrand(double vm, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double Emin = p->P2, Emax = p->P3, Lmin = p->P4, Lmax = p->P5, R = p->P6, psiRz = p->P7, vf = p->P9;
    double E = (psiRz - 0.5 * (pow(vm, 2) + pow(vf, 2))) / psi(0, 0);
//     double L = R * vf * gsl_spline_eval(LcInv, E, lAcc);
    double Rc = Rcirc(E * creal(psi(0, 0)));
    double L = R * vf / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
    if (Emin < E && E < Emax && Lmin < L && L < Lmax) return vm * exp(gsl_spline2d_eval(psdf, E, L, xAcc, yAcc)) - 1;
    else return 0;
}

double PSDF_binned_vm_integral(void * params) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 * p = (struct params_9 *) params;
    
    double Emin = p->P2 * creal(psi(0, 0)), Emax = p->P3 * creal(psi(0, 0)), psiRz = p->P7, vf = p->P9;
    double vm_min2 = (2. * (psiRz - Emax) - pow(vf, 2));
    double vm_max2 = (2. * (psiRz - Emin) - pow(vf, 2));
//     double vm_min2 = 0;
//     double vm_max2 = 2. * psiRz - pow(vf, 2);
    if (vm_max2 <= 0) return 0;
    else {
        if (vm_min2 < 0) vm_min2 = 0;
        gsl_function F;
        F.function = &PSDF_binned_vm_integrand;
        F.params = p;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qags(&F, sqrt(vm_min2), sqrt(vm_max2), 0, tol, nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
    // 	printf("Phi int result: %g\n", result);
        return result;
    }
}

double PSDF_binned_vf_integrand(double vf, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    p->P9 = vf;
    return PSDF_binned_vm_integral(p);
}

double PSDF_binned_vf_integral(void * params) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 * p = (struct params_9 *) params;
    
    double Lmin = p->P4, Lmax = p->P5, R = p->P6, psiRz = p->P7, vRz = p->P8;
    
//     printf("vf = [%g, %g]\n", Lz_min, Lz_max);
    
    gsl_function F;
    F.function = &PSDF_binned_vf_integrand;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
//     gsl_integration_qags(&F, Lz_min / R, Lz_max / R, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_qags(&F, -vRz, vRz, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Phi int result: %g\n", result);
    return result;
}

double PSDF_binned_z_integrand(double z, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = p->P6;
    p->P7 = creal(psi(pow(R, 2), pow(z, 2)));
    p->P8 = sqrt(2. * p->P7);
    return PSDF_binned_vf_integral(p);
}

double PSDF_binned_z_integral(void * params) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 * p = (struct params_9 *) params;
    
    double rmax = p->P1, R = p->P6;
    
    gsl_function F;
    F.function = &PSDF_binned_z_integrand;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, sqrt(pow(rmax, 2) - pow(R, 2)), 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

double PSDF_binned_R_integrand(double R, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    p->P6 = R;
    return R * PSDF_binned_z_integral(p);
}

double PSDF_binned(double Emin, double Emax, double Lmin, double Lmax) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    
    double rmax = 200;
    struct params_9 p = {rmax, Emin, Emax, Lmin, Lmax};
    
    gsl_function F;
    F.function = &PSDF_binned_R_integrand;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, rmax, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return pow(2. * M_PI, 2) * result / 7.5e6;
}



/* ----------------------------------------------------------------------
* Functions used to compute the PSDF reconstructed velocity distributions
* ----------------------------------------------------------------------- */

double pv_num_mag(double v, double R, double psiRz, double rhoR) {
    if (v*v < 2. * psiRz) {
        struct params_9 params = {R, psiRz, v};
        return 4. * M_PI * pow(v, 2) / rhoR * rho_num_c_integral(&params);
    } else return 0;
}

double pv_num_vz_integrand(double vz, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = p->P1;
    double psiRz = p->P2;
    double vr = p->P3;
    double vf = p->P4;
    
    double E = (psiRz - 0.5 * (pow(vr, 2) + pow(vf, 2) + pow(vz, 2))) / psi(0, 0);
    double L = R * vf * gsl_spline_eval(LcInv, E, lAcc);
//     double Rc = Rcirc(E * psi(0, 0));
//     double L = R * vf / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
    
    if (E > 0 && fabs(L) < 1) return exp(gsl_spline2d_eval(psdf, E, L, xAcc, yAcc)) - 1;
    else return 0;
}

double pv_num_vz_integral(void * params) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 * p = (struct params_9 *) params;
    
    double psiRz = p->P2;
    double vr = p->P3;
    double vf = p->P4;
    double vMax = sqrt(2. * psiRz - pow(vr, 2) - pow(vf, 2));
    
    gsl_function F;
    F.function = &pv_num_vz_integrand;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, vMax, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

double pv_num_vr_integrand(double vr, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    p->P3 = vr;
    return pv_num_vz_integral(p);
}

double pv_num_vf(double vf, double R, double psiRz, double rhoR) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    
    if (pow(vf, 2) < 2. * psiRz) {
        double vMax = sqrt(2. * psiRz - pow(vf, 2));
        struct params_9 p = {R, psiRz, 0, vf};
        
        gsl_function F;
        F.function = &pv_num_vr_integrand;
        F.params = &p;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qags(&F, 0, vMax, 0, tol, nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return 4. * result / rhoR;
    } else return 0;
}

double pv_num_vf_integrand(double vf, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    p->P4 = vf;
    return pv_num_vz_integral(p);
}

double pv_num_vr(double vr, double R, double psiRz, double rhoR) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    
    if (pow(vr, 2) < 2. * psiRz) {
        double vMax = sqrt(2. * psiRz - pow(vr, 2));
        struct params_9 p = {R, psiRz, vr, 0};
        
        gsl_function F;
        F.function = &pv_num_vf_integrand;
        F.params = &p;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
        gsl_integration_qags(&F, -vMax, vMax, 0, tol, nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return 2. * result / rhoR;
    } else return 0;
}




/* ---------------------------------
* Parameter initialization function 
* --------------------------------- */

void init_parameters(double* p) {
    A_axi_thin = p[0];
    a_axi_thin = p[1];
    b_axi_thin = p[2];
    A_axi_thick = p[3];
    a_axi_thick = p[4];
    b_axi_thick = p[5];
    A_sph = p[6];
    a_sph = p[7];
    rhoS = pow(10, p[8] + 9);
    rS = p[9];
    gama = p[10];
    alpha = p[11];
    beta = p[12];
    q = p[13];
    vPhi = p[14];
    rPhi2 = p[15];
    
    psiDM0 = 0;
    psiDM0_num = 0;
    if (q == 1) {
        e = 1;
        casin_e = 1;
        sin_e = 1;
    } else {
        e = csqrt(1. - q * q);
        casin_e = casin(e);
        sin_e = creal(casin(e) / e);
    }
    
    if (vPhi > 0) eval_PSDF = &asym_PSDF;
    else {
        vPhi = -vPhi;
        eval_PSDF = &sym_PSDF;
    }
    
    if (q == 1 && alpha == 1 && beta == 3 && gama == 1) {
        if (verbose) printf("Analytic NFW\n");
        rho = &rhoNFW;
        drho_dm2 = &drhoNFW_dm2;
        drho_dz2 = &drhoNFW_dz2;
        d2rho_dz22 = &d2rhoNFW_dz22;
        psiDM = &psiNFW;
        dpsiDM_dR2 = &dpsiNFW_dr2;
        dpsiDM_dz2 = &dpsiNFW_dr2;
        d2psiDM_dz22 = &d2psiNFW_dr22;
        d2psiDM_dz2dR2 = &d2psiNFW_dr22;
    } else if (q == 1 && gama == -1) {
        if (verbose) printf("Analytic BUR\n");
        rho = &rhoBUR;
        drho_dm2 = &drhoBUR_dm2;
        drho_dz2 = &drhoBUR_dz2;
        d2rho_dz22 = &d2rhoBUR_dz22;
        psiDM = &psiBUR;
        dpsiDM_dR2 = &dpsiBUR_dr2;
        dpsiDM_dz2 = &dpsiBUR_dr2;
        d2psiDM_dz22 = &d2psiBUR_dr22;
        d2psiDM_dz2dR2 = &d2psiBUR_dr22;
    } else {
        psiDM = &psiDM_num;
        dpsiDM_dR2 = &dpsiDM_dR2_num;
        dpsiDM_dz2 = &dpsiDM_dz2_num;
        d2psiDM_dz22 = &d2psiDM_dz22_num;
        d2psiDM_dz2dR2 = &d2psiDM_dz2dR2_num;
        if (alpha > 0) {
            if (verbose) printf("Numeric NFW\n");
            rho = &rhoNFW;
            drho_dm2 = &drhoNFW_dm2;
            drho_dz2 = &drhoNFW_dz2;
            d2rho_dz22 = &d2rhoNFW_dz22;
        } else {
            if (verbose) printf("Numeric BUR\n");
            rho = &rhoBUR;
            drho_dm2 = &drhoBUR_dm2;
            drho_dz2 = &drhoBUR_dz2;
            d2rho_dz22 = &d2rhoBUR_dz22;
        }
    }
    
    if (verbose) {
        printf("\t --- Parameter initialization ---\n");
        printf("Thin disc: A = %g, a = %g, b = %g\n", p[0], p[1], p[2]);
        printf("Thick disc: A = %g, a = %g, b = %g\n", p[3], p[4], p[5]);
        printf("Bulge: A = %g, a = %g\n", p[6], p[7]);
        printf("DM: rho_s = %g, r_s = %g, gamma = %g, alpha = %g, beta = %g, q = %g, v_psi = %g, r_a = %g\n", p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]);
        printf("Psi0: %g\n", creal(psi(0,0)));
    }
}



/* --------------------------------------------------------------------
* Functions used to compute the astrophysical direct-detection factors
* -------------------------------------------------------------------- */

double etaV_integrand(double v, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = (p->P1), psiNorm = (p->P2), psiR = (p->P3), ct = (p->P4), st = (p->P5), cf = (p->P6), sf = (p->P7), cT = (p->P8), sT = (p->P9);
    
    /*
    double vx = v * cf * st + Usol + vE * sT * 0.5;
    double vy = v * sf * st + LSR + Vsol + vE * cT;
    double vz = v * ct + Wsol + vE * sT * 0.866;
    */
    double vx = v * cf * st + Usol + vE * (0.9931 * cE - 0.0670 * sE);
    double vy = v * sf * st + LSR + Vsol + vE * (0.1170 * cE + 0.4927 * sE);
    double vz = v * ct + Wsol + vE * (-0.01032 * cE - 0.8676 * sE);
    
    double v2 = pow(vx, 2) + pow(vy, 2) + pow(vz, 2);
    
    double E = (psiR - 0.5 * v2) * psiNorm;
    if (E < 0) return 0;
    double L = vy * R * gsl_spline_eval(LcInv, E, lAcc);
// 	double L = vy * R / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
// 	printf("E = %g\tL = %g\t- psiR: %g, v2: %g, vy: %g, R: %g\n", E, L, psiR, v2, vy, R);
    return pow(v, vPow) * (exp(gsl_spline2d_eval(psdf, E, L, xAcc, yAcc)) - 1);
}

double etaV(void * params) {

    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    
    struct params_9 * p = (struct params_9 *) params;
    double R = (p->P1), psiNorm = (p->P2), psiR = (p->P3), vEsc = (p->P4), vMin = (p->P5), ct = (p->P6), st = (p->P7), cf = (p->P8), sf = (p->P9);
    struct params_9 pp = {R, psiNorm, psiR, ct, st, cf, sf, 0, 0};
    
    double cA = - (cf * st * (Usol + vE * (0.9931 * cE - 0.0670 * sE)) + sf * st * (Vsol + LSR + vE * (0.1170 * cE + 0.4927 * sE)) + ct * (Wsol + vE * (-0.01032 * cE - 0.8676 * sE)));
    double cB = - (pow(Vsol + LSR, 2) + pow(Usol,2 ) + pow(Wsol, 2) + pow(vE, 2) + vE * (sE * (0.9854 * (LSR + Vsol) - 0.1358 * Usol - 1.7352 * Wsol) + cE * (0.234 * (LSR + Vsol) + 1.9862 * Usol - 0.02064 * Wsol - 0.00164952 * vE * sE)));
    double vMax = cA + sqrt(pow(cA, 2) + cB + 2 * psiR);
    
    gsl_function F;
    F.function = &etaV_integrand;
    F.params = &pp;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    if (vMax > vMin) {
        gsl_integration_qags(&F, vMin, vMax, 0, tol, nIntervals, workspace, &result, &abserr);
    } else result = 0;
    gsl_integration_workspace_free(workspace);
    //printf("V int result: %g\n", result);
    return result;
}

double etaPhi_integrand(double f, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    //double cf = cos(f), sf = sin(f);
    double cf = 0, sf = 0;
    sincos(f, &sf, &cf);
    p->P8 = cf;
    p->P9 = sf;
    return etaV(p);
}

double etaPhi(void * params) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 * p = (struct params_9 *) params;

    gsl_function F;
    F.function = &etaPhi_integrand;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, 2. * M_PI, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Phi int result: %g\n", result);
    return result;
}

double etaTheta_integrand(double t, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    //double ct = cos(t), st = sin(t);
    double ct = 0, st = 0;
    sincos(t, &st, &ct);
    p->P6 = ct;
    p->P7 = st;
    double result = st * etaPhi(p);
    return result;
}

double eta(double R, double psiNorm, double psiR, double vEsc, double vMin) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 params = {R, psiNorm, psiR, vEsc, vMin, 0, 0, 0, 0};
    
    gsl_function F;
    F.function = &etaTheta_integrand;
    F.params = &params;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, M_PI, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Theta int result: %g\n", result);
    return result;
}

double *eval_eta(int N, double* results) {
    double psiNorm = 1. / creal(psi(0, 0));
    double psiR = creal(psi(Rsol * Rsol, 0));
    double vEsc = sqrt(2. * psiR);
    double vMax = 1e3;//vEsc + sqrt(pow(LSR + Vsol + vE, 2) + pow(Usol, 2) + pow(Wsol, 2));
    double rhoInv = 1.;// / rho(Rsol * Rsol);
    
    for (int i = 0; i < N; i++) {
        clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
        double vMin = i * vMax / (N - 1.);
        if (vMin == 0) vMin = 1e-3;
        results[3 * i] = vMin;
        vPow = 1;
        results[3 * i + 1] = eta(Rsol, psiNorm, psiR, vEsc, vMin) * rhoInv;
        vPow = 3;
        results[3 * i + 2] = eta(Rsol, psiNorm, psiR, vEsc, vMin) * rhoInv;
        if (verbose) printf("vMin: %g \t-> %g\t%g (computed in %g s)\n", vMin, results[3 * i + 1], results[3 * i + 2], (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
        
    }
    return results;
}

void compute_eta(int v, int n, int m, int N, double* p, double* results) {
    gsl_set_error_handler(&errorHandlerFunc);
    verbose = v;
    init_parameters(p);
    
    interpolate_PSDF(n, m);
    eval_eta(N, results);
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
}

void compute_modulation(int v, int n, int m, int N, int Nmod, double* p, double* results) {
    gsl_set_error_handler(&errorHandlerFunc);
    verbose = v;
    init_parameters(p);
    
    interpolate_PSDF(n, m);
    
    for (int i = 0; i < Nmod; i++) {
        double t = M_PI * i / (Nmod - 1) + 1.15;//1.80394;
        sincos(t, &sE, &cE);
        
        double* r = malloc(3 * N * sizeof(double));
        eval_eta(N, r);
        
        for (int j = 0; j < 3 * N; j++) results[3 * N * i + j] = r[j];
        free(r);
    }
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
}

void compute_pv(int v, int n, int m, int N, double * p, double* results) {
    gsl_set_error_handler(&errorHandlerFunc);
    verbose = v;
    init_parameters(p);
    
    interpolate_PSDF(n, m);
        
    double R2 = Rsol*Rsol, z2 = 0;
    double rhoLoc = rho(R2 + z2 / (q*q));
    double psiRz = psi(R2, z2);
    double vMax = sqrt(2. * psiRz);
    
    for (int i = 0; i < N; i++) {
        double v = 1. * i * vMax / (N - 1);
        struct params_9 params = {Rsol, psiRz, v, 0, 0, 0, 0, 0, 0};
        double Pv = 4. * M_PI * v * v / rhoLoc * rho_num_c_integral(&params);
        if (v < vMax) {
            results[2 * i] = v;
            results[2 * i + 1] = Pv;
        } else {
            results[2 * i] = vMax;
            results[2 * i + 1] = 0;
        }
    }
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
    
}

void compute_pv_components(int v, int n, int m, int N, double * p, double* results) {
    gsl_set_error_handler(&errorHandlerFunc);
    verbose = v;
    init_parameters(p);
    
    interpolate_PSDF(n, m);
        
    double R2 = Rsol*Rsol, z2 = 0;
    double psiRz = psi(R2, z2);
    double vMax = sqrt(2. * psiRz);
    
    results[0] = -vMax;
    results[3 * (N - 1)] = vMax;
    for (int i = 0; i < N; i++) {
        double v = vMax * (2. * i / (N - 1) - 1);
        if (fabs(v) < vMax) {
            results[3 * i] = v;
//             results[3 * i + 1] = pv_num_vf_integral(R2, z2, v);
//             results[3 * i + 2] = pv_num_vr_integral(R2, z2, v);
        }
    }
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
}


void compute_pvs(int v, int n, int m, int N, double * p, double* results_mag, double* results_comp) {
//     gsl_set_error_handler(&errorHandlerFunc);
//     verbose = v;
//     init_parameters(p);
    
//     interpolate_PSDF(n, m);
        
    double R2 = Rsol*Rsol, z2 = 0;
    double rhoLoc = rho(R2 + z2 / (q*q));
    double psiRz = psi(R2, z2);
    double vMax = sqrt(2. * psiRz);
    
    results_mag[2 * (N - 1)] = vMax;
    for (int i = 0; i < N; i++) {
        double v = 1. * i * vMax / (N - 1);
        struct params_9 params = {Rsol, psiRz, v, 0, 0, 0, 0, 0, 0};
        double Pv = 4. * M_PI * v * v / rhoLoc * rho_num_c_integral(&params);
        if (v < vMax) {
            results_mag[2 * i] = v;
            results_mag[2 * i + 1] = Pv;
        }
    }
    
    results_comp[0] = -vMax;
    results_comp[3 * (N - 1)] = vMax;
    for (int i = 0; i < N; i++) {
        double v = vMax * (2. * i / (N - 1) - 1);
        if (fabs(v) < vMax) {
            results_comp[3 * i] = v;
//             results_comp[3 * i + 1] = pv_num_vf_integral(R2, z2, v);
//             results_comp[3 * i + 2] = pv_num_vr_integral(R2, z2, v);
        }
    }
    
//     gsl_spline2d_free(psdf);
//     gsl_spline_free(LcInv);
//     gsl_interp_accel_free(xAcc);
//     gsl_interp_accel_free(yAcc);
//     gsl_interp_accel_free(lAcc);
}


/* ----------------------------------------------------------------------
* Functions used to compute the astrophysical indirect-detection factors
* ---------------------------------------------------------------------- */


double non_enh(double vf, double vm2, double uf, double um2) {
    return 2. * M_PI;
}

double p_wave(double vf, double vm2, double uf, double um2) {
    return 2. * M_PI * (pow(vf - uf, 2) + vm2 + um2);
}

double Sum2_integrand(double um2, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = (p->P1), psiNorm = (p->P2), psiR = (p->P3), vMax = (p->P4), vf = (p->P5), vm2 = (p->P6), uf = (p->P7);
    
    double E = (psiR - 0.5 * (um2 + uf * uf)) * psiNorm;
    double L = R * uf * gsl_spline_eval(LcInv, E, lAcc);
// 	double L = R * uf / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
    
    double result = Sv(vf, vm2, uf, um2) * (exp(gsl_spline2d_eval(psdf, E, L, xAcc, yAcc)) - 1);
    return result;
}

double Sum2(double R, double psiNorm, double psiR, double vMax, double vf, double vm2, double uf) {
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 params = {R, psiNorm, psiR, vMax, vf, vm2, uf, 0, 0};
    
    gsl_function F;
    F.function = &Sum2_integrand;
    F.params = &params;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, 2. * psiR - uf * uf, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

double Suf_integrand(double uf, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = (p->P1), psiNorm = (p->P2), psiR = (p->P3), vMax = (p->P4), vf = (p->P5), vm2 = (p->P6);
    
    double result = Sum2(R, psiNorm, psiR, vMax, vf, vm2, uf);
    return result;
}

double Suf(double R, double psiNorm, double psiR, double vMax, double vf, double vm2) {
// 	clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 params = {R, psiNorm, psiR, vMax, vf, vm2, 0, 0, 0};
    
    gsl_function F;
    F.function = &Suf_integrand;
    F.params = &params;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, -vMax, vMax, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Uf int result: %g in %gs\n", result, (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
    return result;
}

double Svm2_integrand(double vm2, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = (p->P1), psiNorm = (p->P2), psiR = (p->P3), vMax = (p->P4), vf = (p->P5);
    
    double E = (psiR - 0.5 * (vm2 + vf * vf)) * psiNorm;
    double L = R * vf * gsl_spline_eval(LcInv, E, lAcc);
// 	double L = R * vf / (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0)));
    
    double result = Suf(R, psiNorm, psiR, vMax, vf, vm2) * (exp(gsl_spline2d_eval(psdf, E, L, xAcc, yAcc)) - 1);
    return result;
}

double Svm2(double R, double psiNorm, double psiR, double vMax, double vf) {
// 	clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 params = {R, psiNorm, psiR, vMax, vf, 0, 0, 0, 0};
    
    gsl_function F;
    F.function = &Svm2_integrand;
    F.params = &params;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, 2. * psiR - vf * vf, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Vm2 int result: %g in %gs\n", result, (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
    return result;
}


double Svf_integrand(double vf, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double R = (p->P1), psiNorm = (p->P2), psiR = (p->P3), vMax = (p->P4);
    
    double result = Svm2(R, psiNorm, psiR, vMax, vf);
    return result;
}

double Svf(double R, double psiNorm, double psiR, double vMax) {
// 	clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 params = {R, psiNorm, psiR, vMax, 0, 0, 0, 0, 0};
    
    gsl_function F;
    F.function = &Svf_integrand;
    F.params = &params;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, -vMax, vMax, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Vf int result: %g in %gs\n", result, (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
    return result;
}

double Jl1_integrand(double l, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double sf = (p->P1), cf = (p->P2), D = (p->P3), psiNorm = (p->P4);
    double R = D - l * sf;
    double z = l * cf;
    double psiR = creal(psi(R*R, z*z));
    double vMax = sqrt(2. * psiR);
    
    double result = Svf(R, psiNorm, psiR, vMax);
    return result;
}

double Jl2_integrand(double l, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double sf = (p->P1), cf = (p->P2), D = (p->P3), psiNorm = (p->P4);
    double R = l * sf - D;
    double z = l * cf;
    double psiR = creal(psi(R*R, z*z));
    double vMax = sqrt(2. * psiR);
    
    double result = Svf(R, psiNorm, psiR, vMax);
    return result;
}

double Jl(double f, double D, double psiNorm) {
    clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    size_t nIntervals = 1e5;
    double result1, result2, abserr1, abserr2, tol = tolEta;
    double sf, cf;
    sincos(f, &sf, &cf);
    struct params_9 params = {sf, cf, D, psiNorm, 0, 0, 0, 0, 0};
    
    gsl_function F1;
    F1.function = &Jl1_integrand;
    F1.params = &params;
    gsl_integration_workspace *workspace1 = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F1, 0, D - 1e-2, 0, tol, nIntervals, workspace1, &result1, &abserr1);
    gsl_integration_workspace_free(workspace1);
    
    gsl_function F2;
    F2.function = &Jl2_integrand;
    F2.params = &params;
    gsl_integration_workspace *workspace2 = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F2, D + 1e-2, 2. * D, 0, tol, nIntervals, workspace2, &result2, &abserr2);
    gsl_integration_workspace_free(workspace2);
    
    double result = result1 + result2;
    printf("Jl int result: %g in %gs\n", result, (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
    return result;
}

double Jf_integrand(double f, void * params) {
    struct params_9 * p = (struct params_9 *) params;
    double psiNorm = (p->P1), D = (p->P2);
    
    double result = Jl(f, D, psiNorm);
    return result;
}

double Jf(double fMax, double D, double psiNorm, int enhType) {
// 	clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    size_t nIntervals = 1e5;
    double result, abserr, tol = tolEta;
    struct params_9 params = {psiNorm, D, 0, 0, 0, 0, 0, 0, 0};
    
    if (enhType == 0) Sv = non_enh;
    if (enhType == 2) Sv = p_wave;
    
    gsl_function F;
    F.function = &Jf_integrand;
    F.params = &params;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, fMax, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
// 	printf("Jf int result: %g in %gs\n", result, (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
    return result;
}


void compute_PSDF(int v, int n, int m, double* p) {
    gsl_set_error_handler(&errorHandlerFunc);
    verbose = v;
    init_parameters(p);
    
    interpolate_PSDF(n, m);
    
// 	printf("f(E,L): %g\n", exp(gsl_spline2d_eval(psdf, 0.9000000000000000222, 0.1111, xAcc, yAcc)) - 1);
// 	printf("psi0: %g\n", creal(psi(0, 0)));
    
// 	printf("Lc(E): %g\n", gsl_spline_eval(LcInv, 1 - 1e-4, lAcc));
    
// 	double R2 = Rsol, z2 = 0;
// 	double rhoAna = rho(R2 + z2 / (q*q));
// 	double rhoNum = rho_num(R2, z2);
// 	printf("rho: %g, %g, %g\n", rhoAna, rhoNum, rhoNum / rhoAna);
    
    
    
    int N = 20;
    FILE *fp;
    fp = fopen("rho.dat", "w+");
    /**/
    for (int i = 0; i < N + 1; i++) {
        double R2 = (exp(i - N) - exp(- N)) * 1e6 + 0.25;
        fprintf(fp, "%g\t%g\t%g\n", sqrt(R2), creal(rho(R2)), rho_num(R2, 0));
    }
    fclose(fp);
    
    fp = fopen("Lc.dat", "w+");
    N = 1000;
    for (int i = 0; i < N; i++) {
        double E = 1. * i / N;
        double Rc = Rcirc(E * psi(0, 0));
        fprintf(fp, "%g\t%g\t%g\n", E, (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0))), 1. / gsl_spline_eval(LcInv, E, lAcc));
    }
    fclose(fp);
    
    fp = fopen("F.dat", "w+");
    N = 100;
    for (int i = 0; i < N; i++) {
        double E = 1. - pow(10., -4. + 4. * i / (N - 1.));
        double F0 = exp(gsl_spline2d_eval(psdf, E, 0, xAcc, yAcc)) - 1;
        double F1 = exp(gsl_spline2d_eval(psdf, E, 1., xAcc, yAcc)) - 1;
        fprintf(fp, "%g\t%g\t%g\n", E, F0, F1);
    }
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
}

double PSDF_ext(double E, double L, double *p) {
    size_t nIntervals = 1e5;
    double tol = tolPSDF;
    double result, abserr;
    
    verbose = 1;
    gsl_set_error_handler(&errorHandlerFunc);
    
    init_parameters(p);
    psiDM0 = psi(0, 0);
    
    E = E * psiDM0;
    double Rc = Rcirc(E);
    double psiEnv = psi(Rc * Rc, 0);
    L = L * pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0));
    h = 0.05 * psiEnv;
    
    if (verbose) printf("Evaluating PSDF(%g, %g; %g)\n", E, L, psiDM0);
    
    tabulate_z2(E, L, psiEnv, 200);
    
    struct psdf_params params = {E, L, Rc, psiEnv};
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_function F;
    F.params = &params;
    F.function = &PSDF_integrand;
    gsl_integration_qags(&F, 0, M_PI, 0, tol, nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    gsl_spline_free(z2_re);
    gsl_spline_free(z2_im);
    gsl_interp_accel_free(z2_re_acc);
    gsl_interp_accel_free(z2_im_acc);
    
    if (verbose && abserr / result > tol) printf("Warning: error larger then tolerance for PSDF(%g, %g)\n", E, L);
    return result / (2. * pow(M_PI, 2) * sqrt(2));
}


/* ---------------------------------------
* External call for MCMC fit of Vc and Kz
* --------------------------------------- */

double vc2(double R2, double z2) {
    return - 2. * creal(dpsi_dR2(R2, z2));
}

double kz(double R2, double z2) {
    return 2. * csqrt(z2) * creal(dpsi_dz2(R2, z2));
}

double compute_mse(int n_vc, int n_kz, double* vc_data, double* kz_data, double* params) {
    verbose = 0;
    gsl_set_error_handler(&errorHandlerFunc);
    init_parameters(params);
    
    double mse = 0;
    for (int i = 0; i < n_vc; i++) {
        double th = vc2(vc_data[i], 0);
        if (vc_data[i + n_vc] > th) mse += pow((vc_data[i + n_vc] - th), 2) * vc_data[i + 2 * n_vc];
        else mse += pow((vc_data[i + n_vc] - th), 2) * vc_data[i + 3 * n_vc];
    }
    
    for (int i = 0; i < n_kz; i++) {
        double th = kz(kz_data[i], kz_data[i + n_kz]);
        if (kz_data[i + 2 * n_kz] > th) mse += pow((kz_data[i + 2 * n_kz] - th), 2) * kz_data[i + 3 * n_kz];
        else mse += pow((kz_data[i + 2 * n_kz] - th), 2) * kz_data[i + 4 * n_kz];
    }
    
    return 0.5 * mse;
}

double* compute_psi(int Nh, int Nv, double x) {
    verbose = 0;
    gsl_set_error_handler(&errorHandlerFunc);
    
    double* results = malloc((Nh + Nv) * sizeof(double));
    
    double a = log10(0.1);
    double b = (log10(100) - a) / (Nh - 1);
    for (int i = 0; i < Nh; i++) {
        double R = pow(10, b * i + a);
        results[i] = psi_exp_num(pow(R, 2), 1e-3 * x, 1, x, 0);
        printf("psi(%g, %g) = %g\n", R, 0., results[i]);
    }
    
    a = 0.01;
    b = (10 - a) / (Nv - 1);
    for (int i = 0; i < Nv; i++) {
        double z = b * i + a;
        results[i + Nh] = psi_exp_num(pow(2.5, 2), z * x, 1, x, 0);
        printf("psi(%g, %g) = %g\n", 2.5, z, results[i + Nh]);
    }
    
    return results;
}

double* compute_moments(int ver, int n, int m, int N, int Nr, double * Rpts, double * p) {
    verbose = ver;
    gsl_set_error_handler(&errorHandlerFunc);
    init_parameters(p);	
    interpolate_PSDF(n, m);
    
    double vMax = sqrt(2. * psi(0, 0));
    double* results = malloc((2 + 3 * Nr) * N * sizeof(double));
    for (int i = 0; i < N * (2 + 3 * Nr); i++) results[i] = 0;
    for (int i = 0; i < N; i++) results[i] = vMax * i / (N - 1.);
    for (int i = 0; i < N; i++) results[i + N * (1 + Nr)] = vMax * (2. * i / (N - 1) - 1);
    
    for (int x = 0; x < Nr; x++) {
        double R2 = pow(Rpts[x], 2), z2 = 0;
        double rhoLoc = rho(R2 + z2 / (q*q));
        double psiRz = psi(R2, z2);
    
        for (int i = 0; i < N; i++) {
            double v = results[i];
            struct params_9 params = {Rsol, psiRz, v, 0, 0, 0, 0, 0, 0};
            if (v < vMax) results[N * (1 + x) + i] = 4. * M_PI * v * v / rhoLoc * rho_num_c_integral(&params);
        }
        
        for (int i = 0; i < N; i++) {
            double v = results[i + N * (1 + Nr)];
            if (fabs(v) < vMax) {
//                 results[N * (2 + Nr + 2 * x) + i] = pv_num_vf_integral(R2, z2, v);
//                 results[N * (3 + Nr + 2 * x) + i] = pv_num_vr_integral(R2, z2, v);
            }
        }
    }
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
    
    return results;
}

void analysis(int ver, int n, int m, int N, int Nr, double * Rpts, double * p) {
    verbose = ver;
    gsl_set_error_handler(&errorHandlerFunc);
    init_parameters(p);	
    interpolate_PSDF(n, m);
    
    double vMax = sqrt(2. * psi(0, 0));
    double* results = malloc((2 + 3 * Nr) * N * sizeof(double));
    for (int i = 0; i < N * (2 + 3 * Nr); i++) results[i] = 0;
    for (int i = 0; i < N; i++) results[i] = vMax * i / (N - 1.);
    for (int i = 0; i < N; i++) results[i + N * (1 + Nr)] = vMax * (2. * i / (N - 1) - 1);
    
    for (int x = 0; x < Nr; x++) {
        double R2 = pow(Rpts[x], 2), z2 = 0;
        double rhoLoc = rho(R2 + z2 / (q*q));
        double psiRz = psi(R2, z2);
    
        for (int i = 0; i < N; i++) {
            double v = results[i];
            struct params_9 params = {Rsol, psiRz, v, 0, 0, 0, 0, 0, 0};
            if (v < vMax) results[N * (1 + x) + i] = 4. * M_PI * v * v / rhoLoc * rho_num_c_integral(&params);
        }
        
        for (int i = 0; i < N; i++) {
            double v = results[i + N * (1 + Nr)];
            if (fabs(v) < vMax) {
//                 results[N * (2 + Nr + 2 * x) + i] = pv_num_vf_integral(R2, z2, v);
//                 results[N * (3 + Nr + 2 * x) + i] = pv_num_vr_integral(R2, z2, v);
            }
        }
    }
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
}

void print_integrand(char const *fileName, double complex R2, double complex z2, double Rd, double zd, int derivative) {
    FILE *fp;
    int n = 1000;
    
    switch(derivative) {
        case 0:
            psi_exp_num_kernel = &psi_exp_kernel_0;
            break;
        case 1:
            psi_exp_num_kernel = &psi_exp_kernel_1;
            break;
        case 2:
            psi_exp_num_kernel = &psi_exp_kernel_2;
            break;
        case 3:
            psi_exp_num_kernel = &psi_exp_kernel_3;
            break;
        case 4:
            psi_exp_num_kernel = &psi_exp_kernel_4;
            break;
    }
    
    struct psi_exp_params params = {csqrt(R2), csqrt(z2), 1. / Rd, 1. / zd, 13.5686};
    fp = fopen(fileName, "w+");
    for (int i = 0; i < n; i++) {
        double R = (exp((i - n) / 100.) - exp(- n / 100.)) * exp(1 / 100.) * 1e2 + 1e-8;
        fprintf(fp, "%g\t%g\n", R, psi_exp_integrand1_re(R, &params));
    }
    fclose(fp);
}

void print_Pv(char const *fileName, int n, int m, int N, double * p) {
    gsl_set_error_handler(&errorHandlerFunc);
    init_parameters(p);	
    interpolate_PSDF(n, m);
        
    double R2 = 0, z2 = pow(0.3 * rS, 2);
    double rhoAna = rho(R2 + z2 / (q*q));
    double rhoNum = rho_num(R2, z2);
    printf("rho: %g, %g, %g\n", rhoAna, rhoNum, rhoNum / rhoAna);
    
    
    double rhoLoc = rho(R2 + z2 / (q*q));
    double psiRz = psi(R2, z2);
    double vMax = sqrt(2. * psiRz);
    
    /*
    */
    FILE *fp;
    fp = fopen(fileName, "w+");
    for (int i = 0; i < n; i++) {
        double v = 1. * i * vMax / n;
        struct params_9 params = {sqrt(R2), psiRz, v, 0, 0, 0, 0, 0, 0};
        double Pv = 4. * M_PI * v * v / rhoLoc * rho_num_c_integral(&params);
        fprintf(fp, "%g\t%g\n", v, Pv);
    }
    fprintf(fp, "%g\t%g\n", vMax, 0.);
    fclose(fp);
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
    
}

void print_rho(char const *fileName, int n) {
    FILE *fp;
    fp = fopen(fileName, "w+");
    
    for (int i = 0; i < n; i++) {
        double R = (exp((i - n) / 100.) - exp(- n / 100.)) * exp(1 / 100.) * 1e2 + 1e-8;
        fprintf(fp, "%g\t%g\t%g\n", R, creal(rho(R*R)), rho_num(R, 0));
    }
    fclose(fp);
}

void * PSDF_binned_thread(void *vargp) { 
    // Store the value argument passed to this thread 
    double *args = (double *) vargp; 
    double *r = malloc(sizeof(double)); 
    *r = PSDF_binned(args[0], args[1], args[2], args[3]);
    if (verbose) printf("E = [%g, %g], L = [%g, %g] -> %g\n", args[0], args[1], args[2], args[3], *r);
    return (void *) r;
} 

void HQ_model(int v, int n, int m, int Nv, int Nr, double * rpts, double * p, char const *name) {
    verbose = v;
    gsl_set_error_handler(&errorHandlerFunc);
    init_parameters(p);
    interpolate_PSDF(n, m);
    
    char folder[128], name_vel_mag[128], name_vel_R[128], name_vel_f[128], name_psdf[128], name_lz[128], name_rho[128];
    sprintf(name_vel_mag, "%svelocity_mag.dat", name);
    sprintf(name_vel_R, "%svelocity_R.dat", name);
    sprintf(name_vel_f, "%svelocity_f.dat", name);
    sprintf(name_psdf, "%sdf.dat", name);
    sprintf(name_lz, "%slz.dat", name);
    sprintf(name_rho, "%srho.dat", name);
    
    double z2 = 0;
    double vMax = sqrt(2. * creal(psi(0,0)));
    
    double vel_mag[Nv * (1 + Nr)];
    double vel_R[Nv * (1 + Nr)];
    double vel_f[Nv * (1 + Nr)];
    
    for (int i = 0; i < Nv; i++) {
        vel_mag[i * (1 + Nr)] = vMax * i / (Nv - 1);
        vel_R[i * (1 + Nr)] = vMax * (2. * i / (Nv - 1) - 1);
        vel_f[i * (1 + Nr)] = vMax * (2. * i / (Nv - 1) - 1);
    }
    
    
    
    if (verbose) printf("Computing velocities...\n");
    for (int j = 0; j < Nr; j++) {
        double R2 = pow(rpts[j], 2);
        double psiRz = creal(psi(R2, z2));
        double rhoR = rho(R2 + z2 / (q*q));
        for (int i = 0; i < Nv; i++) {
//             if (verbose) printf("Computing velocities: %d (v_mag = %g, v_R/f = %g)\n", i, vel_mag[i * (1 + Nr)], vel_R[i * (1 + Nr)]);
            vel_mag[i * (1 + Nr) + j + 1] = pv_num_mag(vel_mag[i * (1 + Nr)], rpts[j], psiRz, rhoR);
            vel_R[i * (1 + Nr) + j + 1] = pv_num_vr(vel_R[i * (1 + Nr)], rpts[j], psiRz, rhoR);
            vel_f[i * (1 + Nr) + j + 1] = pv_num_vf(vel_f[i * (1 + Nr)], rpts[j], psiRz, rhoR);
        }
    }
    
    FILE *f_vel_mag, *f_vel_R, *f_vel_f;
    f_vel_mag = fopen(name_vel_mag, "w+");
    f_vel_R = fopen(name_vel_R, "w+");
    f_vel_f = fopen(name_vel_f, "w+");
    for (int i = 0; i < Nv; i++) {
        for (int j = 0; j < Nr + 1; j++) {
            fprintf(f_vel_mag, "%g\t", vel_mag[i * (1 + Nr) + j]);
            fprintf(f_vel_R, "%g\t", vel_R[i * (1 + Nr) + j]);
            fprintf(f_vel_f, "%g\t", vel_f[i * (1 + Nr) + j]);
        }
        fprintf(f_vel_mag, "\n");
        fprintf(f_vel_R, "\n");
        fprintf(f_vel_f, "\n");
    }
    fclose(f_vel_mag);
    fclose(f_vel_R);
    fclose(f_vel_f);
    
    
    
    if (verbose) printf("Computing binned PSDF...\n");
    /*FILE *f_psdf;
    f_psdf = fopen(name_psdf, "w+");
    int N_psdf = 20, M_psdf = 5;
    double E_pts[N_psdf];
    double Lz_pts[2 * M_psdf - 1];
    for (int i = 0; i < N_psdf; i++)  E_pts[i] = pow(2 - 1e-4, 1. * i / (N_psdf - 1)) - 1;
    for (int j = 0; j < 2 * M_psdf - 1; j++) Lz_pts[j] = 2. * j / (2. * M_psdf - 2.) - 1.;
    for (int i = 0; i < N_psdf - 1; i++) {
        for (int j = 0; j < 2 * M_psdf - 2; j++) {
            printf("E = [%g, %g], L = [%g, %g]", E_pts[i], E_pts[i + 1], Lz_pts[j], Lz_pts[j + 1]);
            double df_bin = PSDF_binned(E_pts[i], E_pts[i + 1], Lz_pts[j], Lz_pts[j + 1]);
//             double args[] = {E_pts[i], E_pts[i + 1], Lz_pts[j], Lz_pts[j + 1]};
            
            printf(" -> %g\n", df_bin);
            fprintf(f_psdf, "%g\t", df_bin);
        } fprintf(f_psdf, "\n");
    } fclose(f_psdf);
    */
    /*
    int N_psdf = 20, M_psdf = 5;
    double E_pts[N_psdf], Lz_pts[2 * M_psdf - 1];
    for (int i = 0; i < N_psdf; i++)  E_pts[i] = pow(2 - 1e-4, 1. * i / (N_psdf - 1)) - 1;
    for (int j = 0; j < 2 * M_psdf - 1; j++) Lz_pts[j] = 2. * j / (2. * M_psdf - 2.) - 1.;
    
    int N_threads = 2 * (N_psdf - 1) * (M_psdf - 1);
    pthread_t thread[N_threads];
    for (int i = 0; i < N_threads; i++) {
//         double args[] = {E_pts[i / (2 * M_psdf - 2)], E_pts[i / (2 * M_psdf - 2) + 1], Lz_pts[i % (2 * M_psdf - 2)], Lz_pts[i % (2 * M_psdf - 2) + 1]};
        double *args = malloc(4 * sizeof(double));
        args[0] = E_pts[i / (2 * M_psdf - 2)];
        args[1] = E_pts[i / (2 * M_psdf - 2) + 1];
        args[2] = Lz_pts[i % (2 * M_psdf - 2)];
        args[3] = Lz_pts[i % (2 * M_psdf - 2) + 1];
//         printf("%d (%d, %d): E = [%g, %g], Lz = [%g, %g]\n", i, i / (2 * M_psdf - 2), i % (2 * M_psdf - 2), args[0], args[1], args[2], args[3]);
        pthread_create(&thread[i], NULL, PSDF_binned_thread, (void*)args);
    }
    
    FILE *f_psdf;
    f_psdf = fopen(name_psdf, "w+");
    for(int i = 0; i < N_threads; i++) {
        void *result;
        pthread_join(thread[i], (void **) &result);
        fprintf(f_psdf, "%g\t", * (double *) result);
    } fclose(f_psdf);
    */
    
    
    if (verbose) printf("Computing Lz...\n");
    FILE *f_lz;
    f_lz = fopen(name_lz, "w+");
    for (int i = 0; i < Nv; i++) {
        double E = pow(2 - 1e-4, 1. * i / (Nv - 1)) - 1;//1. * i / N;
        double Rc = Rcirc(E * psi(0, 0));
        fprintf(f_lz, "%g\t%g\t%g\t%g\n", E, E * creal(psi(0, 0)), (pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0))), 1. / gsl_spline_eval(LcInv, E, lAcc));
    } fclose(f_lz);
    
    
    
    if (verbose) printf("Computing rho...\n");
    FILE *fp;
    fp = fopen(name_rho, "w+");
    for (int i = 0; i < Nv + 1; i++) {
        double R2 = (exp(i - Nv) - exp(- Nv)) * 1e6 + 0.25;
        fprintf(fp, "%g\t%g\t%g\n", sqrt(R2), creal(rho(R2)), rho_num(R2, 0));
    } fclose(fp);
    
    gsl_spline2d_free(psdf);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(xAcc);
    gsl_interp_accel_free(yAcc);
    gsl_interp_accel_free(lAcc);
    
}

void * myThreadFun(void *vargp) { 
    // Store the value argument passed to this thread 
    double * val = (double *) vargp; 
    double *r = malloc(sizeof(double)); 
    *r = sqrt(val[1]);
//     printf("%g -> %g\n", *val, *r); 
    return (void *) r;
} 

/* -----------------------------------------------------------------------------------------
* Main function that I was using for testing the code (some lines might still be useful...)
* ----------------------------------------------------------------------------------------- */

int main (int argc, char **argv) {
    verbose = 1;
    gsl_set_error_handler(&errorHandlerFunc);
    clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
    
    double params[16] = {5e10, 3.6, 0.5, 0, 1., 1., 1e11, 1., -2., 13., 1., 1., 3., 1., 0, 16.};
    compute_PSDF(1, 100, 20, params);
    /*
    init_parameters(params);
    
    double complex R2 = 7.3 - 2.7*I, z2 = 100.3 + 5.7*I;
    double complex dpsi = dpsi_dz2(R2, z2);
    double complex r = d2rho_dz22(R2, z2) * cpow(dpsi, -2) - drho_dz2(R2, z2) * d2psi_dz22(R2, z2) * cpow(dpsi, -3);
    printf("R2 = (%g, %gI), z2 = (%g, %gI): (%g, %gI)\n", creal(R2), cimag(R2), creal(z2), cimag(z2), creal(r), cimag(r));
    
    
    double E = 0.9 * psi(0, 0);
    double Rc = Rcirc(E);
    double psiEnv = psi(Rc * Rc, 0);
    double L = 0.1 * pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0));
    h = 0.05 * psiEnv;
    struct psdf_params p = {E, L, Rc, psiEnv};
    printf("df(E=%g, Lz=%g) = %g\n", E, L, PSDF_integrand(M_PI - 1e-2, &p));
    
    
    double E = 0.999 * creal(psi(0, 0));
    double Rc = Rcirc(E);
    printf("Rc: %g\n", Rc);
    double L = 1. * pow(Rc, 2) * sqrt(-2. * dpsi_dR2(Rc * Rc, 0));
    printf("E: %.16g, Lz: %.16g\n", E, L);
    //double complex R2 = 1.18+1.3*I;
    //double complex xi = pow(L, 2) / (2. * R2) + E;
    double complex xi = 573185.0808622467 - 6.2906290463055e-05*I;
    printf("xi: %.16g, %.16g)\n", creal(xi), cimag(xi));
    double complex z2 = psi_inverse(xi, E, L, 1e3, 1e3);
    printf("z2 = (%g, %g)\n", creal(z2), cimag(z2));
    */
    printf("Done in %g s\n", (clock() / (CLOCKS_PER_SEC / 1000) - start) / 1000.);
    
    return 0;
}
