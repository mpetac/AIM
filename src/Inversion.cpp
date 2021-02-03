#include "Inversion.hpp"

/**
 * @param model Galactic model
 * @param N_E Number of relative energy interpolation points
 * @param N_Lz Number of angular momentum interpolation points
 * @param N_Lc Number of maximum angular momentum interpolation points
 * @param tolerance_F Tolerance used for computing the PSDF
 * @param verobse Verbose output
 */

Inversion::Inversion(Model *model, int N_E, int N_Lz, double tolerance_F, bool verbose) {
    Inversion::model = model;
    Inversion::tolerance_F = tolerance_F;
    Inversion::verbose = verbose;
    if (verbose) gsl_set_error_handler(&Inversion::GSL_error_func);
    else gsl_set_error_handler(&Inversion::GSL_error_func_silent);
    
    time_t tStart = time(NULL);
    
    double Epts[N_E];
    double Lzpts[2 * N_Lz - 1];
    double Fpts[N_E * (2 * N_Lz - 1)];
    
    double base = 2.;
    double logEmin = std::log10(2.) / std::log10(base), logEmax = std::log10(1. + 1e-4) / std::log10(base);
    for (int i = 0; i < N_E; i++) {
        Epts[i] = 2. - std::pow(base, logEmin + (logEmax - logEmin) * i / (N_E - 1));
    }
    
    for (int i = 0; i < 2 * N_Lz - 1; i++) {
        Lzpts[i] = 2. * i / (2. * N_Lz - 2.) - 1;
    }
    
    Inversion::tabulate_F(N_E, N_Lz, Epts, Lzpts, Fpts);    
    
    Inversion::EAcc = gsl_interp_accel_alloc();
    Inversion::LzAcc = gsl_interp_accel_alloc();
    Inversion::F = gsl_spline2d_alloc(gsl_interp2d_bilinear, N_E, 2 * N_Lz - 1);
    gsl_spline2d_init(Inversion::F, Epts, Lzpts, Fpts, N_E, 2 * N_Lz - 1);
    
    double dt = difftime(time(NULL), tStart);
    if (Inversion::verbose) std::cout << "Inversion computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
}

Inversion::~Inversion() {
    gsl_spline2d_free(F);
    gsl_interp_accel_free(EAcc);
    gsl_interp_accel_free(LzAcc);
}

/**
 * @param N_E Number of interpolation points for relative energy
 * @param N_Lz Number of interpolation points for angular momentum
 * @param Epts Array of interpolation points for the relative energy
 * @param Lzpts Array of interpolation points for the angular momentum
 * @param Fpts Array where the tabulated PSDF is stored
 */

void Inversion::tabulate_F(int N_E, int N_Lz, double* Epts, double* Lzpts, double* Fpts) {
    if(Inversion::model->psi(1., 0, 1.) == Inversion::model->psi(0, 1., 1.)) {
        if (Inversion::verbose) std::cout << "Computing Eddington's inversion..." << std::endl;
        Inversion::h = 0;
        std::vector<std::future<double>> vals(N_E);
        for (int i = 0; i < N_E; i++) {
            //double *params = new double[2];
            //params[0] = Epts[i];
            //params[1] = 0;
            //vals[i] = std::async(std::launch::async, &Inversion::F_even, this, params);
            vals[i] = std::async(std::launch::async, &Inversion::F_eddington, this, Epts[i]);
        }
        for (int i = 0; i < N_E; i++) {
            double val = std::log(1. + vals[i].get());
            if (val < 0) {
                if (Inversion::verbose) std::cout << "PSDF-even negative: " << Epts[i] << " -> " << val << std::endl;
                val = 0;
            }
            //std::cout << "PSDF computed: " << Epts[i] << " -> " << val << std::endl;
            for (int j = 0; j < N_Lz; j++) {
                Fpts[(N_Lz - 1 + j) * N_E + i] = val;
                Fpts[(N_Lz - 1 - j) * N_E + i] = val;
            }
        }
    } else {
        if (Inversion::verbose) std::cout << "Computing axisymmetric inversion..." << std::endl;
        std::vector<std::future<double>> vals_even(N_E * N_Lz);
        std::vector<std::future<double>> vals_odd(N_E * N_Lz);
        for (int i = 0; i < N_E; i++) {
            for (int j = 0; j < N_Lz; j++) {
                double *params = new double[2];
                params[0] = Epts[i];
                params[1] = Lzpts[j + N_Lz - 1];
                // for single threaded execution change to "std::launch::deferred"
                vals_even[i * N_Lz + j] = std::async(std::launch::async, &Inversion::F_even, this, params);
                vals_odd[i * N_Lz + j] = std::async(std::launch::async, &Inversion::F_odd, this, params);
            }
        }
        
        for (int i = 0; i < N_E; i++) {
            for (int j = 0; j < N_Lz; j++) {
                double val_even = vals_even[i * N_Lz + j].get();
                double val_odd = vals_odd[i * N_Lz + j].get();
                if (val_even < 0) {
                    if (Inversion::verbose) std::cout << "Warning! f+ negative: " << Epts[i] << ", " << Lzpts[j + N_Lz - 1] << " -> " << val_even << ", " << val_odd << std::endl;
                    val_even = 0;
                }
                if (std::abs(val_odd) > val_even) {
                    if (Inversion::verbose) std::cout << "Warning! f- larger then f+: " << Epts[i] << ", " << Lzpts[j + N_Lz - 1] << " -> " << val_even << ", " << val_odd << std::endl;
                    val_odd = val_even * (1. - 2. * int(std::signbit(val_odd)));
                }
//                 std::cout << "PSDF computed: " << Epts[i] << ", " << Lzpts[j + N_Lz - 1] << " -> " << val_even << ", " << val_odd << std::endl;
                Fpts[(N_Lz - 1 + j) * N_E + i] = std::log(1. + val_even + val_odd);
                Fpts[(N_Lz - 1 - j) * N_E + i] = std::log(1. + val_even - val_odd);
            }
        }
    }
}

/**
 * @param t Angle parametrizing the position along the contour
 * @param params Struct with data needed for evaluating the contour integral
 */

double F_even_integrand(double t, void *params) {
    inversion_params *p = (inversion_params *) params;
    Model *model = (Model *) p->model;
    InversionInterp *z2interp = (InversionInterp *) p->z2interp;
    
    std::complex<double> xi = 0.5 * p->psiEnv * (1. + std::cos(t) + 2.i * p->h * std::sin(t));
    if (std::abs(xi) < model->psi0 * 1e-5) return 0;
    std::complex<double> dxi = 0.5 * p->psiEnv * (1.i * std::sin(t) + 2. * p->h * std::cos(t));
    std::complex<double> jac = std::pow(xi - p->E, -0.5);
    std::complex<double> R2 = std::pow(p->Lz, 2) / (2. * (xi - p->E));
    std::complex<double> z2 = model->psi_inverse(xi, p->E, p->Lz, z2interp->z2_eval(t));
    std::complex<double> drho_dpsi = model->rho_d2psi2(R2, z2, std::sqrt(R2 + z2));
    
    //std::cout << R2 << ", " << z2 << std::endl;
    //std::cout << dxi << ", " << jac << ", " << drho_dpsi << std::endl;
    
    return std::real(dxi * jac * drho_dpsi);
}

/**
 * @param params Values of the relative energy and angular momentum (expressed as a vector in E-Lz plane)
 */

double Inversion::F_even(double* params) {
//     if (Inversion::verbose) std::cout << "Computing F-even: " << params[0] << ", " << params[1] << std::endl;
    
    if (params[0] == 0) return 0;
    
    double result, abserr;
    double E = params[0] * Inversion::model->psi0;
    double Rc = Inversion::model->Rcirc(E);
    double Rc2 = std::pow(Rc, 2);
    double psiEnv = std::real(Inversion::model->psi(Rc2, 0, Rc));
    double Lz = params[1] * std::pow(Rc, 2) * std::sqrt(-2. * std::real(Inversion::model->psi_dR2(Rc2, 0, Rc)));
    
    InversionInterp z2interp(Inversion::model, E, Lz, psiEnv, Inversion::h, Inversion::nInterp);
    inversion_params p = {Inversion::model, &z2interp, E, Lz, Rc, psiEnv, Inversion::h};
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Inversion::nIntervals);
    gsl_function F;
    F.params = &p;
    F.function = &F_even_integrand;
    gsl_integration_qags(&F, 0, M_PI, 0, Inversion::tolerance_F, Inversion::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return result * Inversion::result_fact_even;
}

/**
 * @param t Angle parametrizing the position along the contour
 * @param params Struct with data needed for evaluating the contour integral
 */

double F_odd_integrand(double t, void *params) {
    inversion_params *p = (inversion_params *) params;
    Model *model = (Model *) p->model;
    InversionInterp *z2interp = (InversionInterp *) p->z2interp;
    
    std::complex<double> xi = 0.5 * p->psiEnv * (1. + std::cos(t) + 2.i * p->h * std::sin(t));
    if (std::abs(xi) < model->psi0 * 1e-5) return 0;
    std::complex<double> dxi = 0.5 * p->psiEnv * (1.i * std::sin(t) + 2. * p->h * std::cos(t));
    std::complex<double> jac = std::pow(xi - p->E, -0.5);
    std::complex<double> R2 = std::pow(p->Lz, 2) / (2. * (xi - p->E));
    std::complex<double> z2 = model->psi_inverse(xi, p->E, p->Lz, z2interp->z2_eval(t));
    std::complex<double> drhovphi_dpsi = model->rho_vphi_d2psi2(R2, z2, std::sqrt(R2 + z2));
    return std::real(dxi * jac * drhovphi_dpsi);
}

/**
 * @param params Values of the relative energy and angular momentum (expressed as a vector in E-Lz plane)
 */

double Inversion::F_odd(double* params) {
//     if (Inversion::verbose) std::cout << "Computing F-odd: " << params[0] << ", " << params[1] << std::endl;
    
    if (!Inversion::model->is_rotating() || params[0] == 0 || params[1] == 0) return 0;
    
    double result, abserr;
    double E = params[0] * Inversion::model->psi0;
    double Rc = Inversion::model->Rcirc(E);
    double Rc2 = std::pow(Rc, 2);
    double psiEnv = std::real(Inversion::model->psi(Rc2, 0, Rc));
    double Lz = params[1] * std::pow(Rc, 2) * std::sqrt(-2. * std::real(Inversion::model->psi_dR2(Rc2, 0, Rc)));
    
    InversionInterp z2interp(Inversion::model, E, Lz, psiEnv, Inversion::h, Inversion::nInterp);
    inversion_params p = {Inversion::model, &z2interp, E, Lz, Rc, psiEnv, Inversion::h};
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Inversion::nIntervals);
    gsl_function F;
    F.params = &p;
    F.function = &F_odd_integrand;
    gsl_integration_qags(&F, 0, M_PI, 0, Inversion::tolerance_F, Inversion::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return result * Inversion::result_fact_odd;
}


/**
 * @param r Radial distance
 * @param params Struct with data needed for evaluating the contour integral
 */

double F_eddington_integrand(double r, void *params) {
    inversion_eddington_params *p = (inversion_eddington_params *) params;
    Model *model = (Model *) p->model;
    
    double r2 = std::pow(r, 2);
    double psi_r = std::real(model->psi(r2, 0, r));
    double dpsi_dr2 = std::real(model->psi_dR2(r2, 0, r));
    
    double jac = std::pow(p->E - psi_r, -0.5);
    double d2rho_dpsi2 = std::real(model->rho_d2psi2(r2, 0, r));
    
    return -2. * r * dpsi_dr2 * d2rho_dpsi2 * jac;
}

/**
 * @param E Values of the relative energy
 */

double Inversion::F_eddington(double E) {
    if (E == 0) return 0;
    
    double result, abserr;
    double rE = std::sqrt(std::real(Inversion::model->psi_inverse(E * Inversion::model->psi0, 0, 0, 1)));
    inversion_eddington_params p = {Inversion::model, std::real(Inversion::model->psi(std::pow(rE, 2), 0, rE))};
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Inversion::nIntervals);
    gsl_function F;
    F.params = &p;
    F.function = &F_eddington_integrand;
    gsl_integration_qagiu(&F, rE, 0, Inversion::tolerance_F, Inversion::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return result * Inversion::result_fact_even;
}


/**
 * @param E Value of the relative energy
 * @param Lz Value of the angular momentum
 */

double Inversion::eval_F(double E, double Lz) {
    if (E < 0 || E > 1 || Lz < -1 || Lz > 1) {
        if (Inversion::verbose) std::cout << "Warning! PSDF eval out of range: " << E << ", " << Lz << std::endl;
        return 0;
    }
    double f = std::exp(gsl_spline2d_eval(F, E, Lz, EAcc, LzAcc)) - 1.;
    return f;
}

/**
 * @param reason Reason for raising the error
 * @param file Name of the file in which the error occured
 * @param line Line in which the error occured
 * @param gsl_errno GSL error code
 */

void Inversion::GSL_error_func(const char* reason, const char* file, int line, int gsl_errno) {
    std::cout << " -> GSL error #" << gsl_errno << " in line " << line << " of " << file <<": " << gsl_strerror(gsl_errno) << std::endl;
}

/**
 * @param reason Reason for raising the error
 * @param file Name of the file in which the error occured
 * @param line Line in which the error occured
 * @param gsl_errno GSL error code
 */

void Inversion::GSL_error_func_silent(const char* reason, const char* file, int line, int gsl_errno) {
}

