#include "Inversion.hpp"

Inversion::Inversion(Model *model, int N_E, int N_Lz, int N_Lc, double tolerance_F, bool verbose) {
    Inversion::model = model;
    Inversion::tolerance_F = tolerance_F;
    Inversion::verbose = verbose;
    Inversion::psi0 = std::real(Inversion::model->psi(0, 0, 0));
    
    gsl_set_error_handler(&GSL_error_func);
    
    double Epts[N_E];
    double Lzpts[N_Lz];
    double Fpts[N_E * (2 * N_Lz - 1)];
    
    for (int i = 0; i < N_E; i++) {
        Epts[i] = std::pow(2. - 1e-4, 1. * i / (N_E - 1)) - 1;
    }
    
    for (int i = 0; i < N_Lz; i++) {
        Lzpts[i] = 2. * i / (2. * N_Lz - 2.) - 1;
    }
    
    Inversion::tabulate_F(N_E, N_Lz, Epts, Lzpts, Fpts);    
    
    Inversion::EAcc = gsl_interp_accel_alloc();
    Inversion::LzAcc = gsl_interp_accel_alloc();
    Inversion::F = gsl_spline2d_alloc(gsl_interp2d_bilinear, N_E, 2 * N_Lz - 1);
    gsl_spline2d_init(Inversion::F, Epts, Lzpts, Fpts, N_E, 2 * N_Lz - 1);
    
    
    
    double Lcpts[N_Lc];
    double LcI[N_Lc];
    
    double epsilon = 1.02;
    for (int i = 0; i < N_Lc; i++) {
        Lcpts[N_Lc - 1 - i] = 1. - std::pow(epsilon, i) / std::pow(epsilon, N_Lc - 1);
    }
    
    for (int i = 0; i < N_Lc; i++) {
        double Rc = model->Rcirc(Lcpts[i] * Inversion::psi0);
        double psi_dR2_Rc = std::real(model->psi_dR2(std::pow(Rc, 2), 0, 0));
        LcI[i] = 1. / (std::pow(Rc, 2) * std::sqrt(-2. * psi_dR2_Rc));
    }
    
    Inversion::LcAcc = gsl_interp_accel_alloc();
    Inversion::LcInv = gsl_spline_alloc(gsl_interp_cspline, N_Lc);
    gsl_spline_init(Inversion::LcInv, Lcpts, LcI, N_Lc);
}

Inversion::~Inversion() {
    gsl_spline2d_free(F);
    gsl_spline_free(LcInv);
    gsl_interp_accel_free(EAcc);
    gsl_interp_accel_free(LzAcc);
    gsl_interp_accel_free(LcAcc);
}


void Inversion::tabulate_F(int N_E, int N_Lz, double* Epts, double* Lzpts, double* Fpts) {
    if(Inversion::model->psi(1., 0, 1.) == Inversion::model->psi(0, 1., 1.)) {
        std::vector<std::future<double>> vals(N_E);
        for (int i = 0; i < N_E; i++) {
            double params[2] = {Epts[i], 0};
            vals[i] = std::async(&Inversion::F_even, this, params);
        }
        for (int i = 0; i < N_E; i++) {
            double val = std::log(1. + vals[i].get());
            for (int j = 0; j < N_Lz; j++) {
                Fpts[(N_Lz - 1 + j) * N_E + i] = val;
                Fpts[(N_Lz - 1 - j) * N_E + i] = val;
            }
        }
    } else {
        std::vector<std::future<double>> vals_even(N_E * N_Lz);
        std::vector<std::future<double>> vals_odd(N_E * N_Lz);
        for (int i = 0; i < N_E; i++) {
            for (int j = 0; j < N_Lz; j++) {
                double params[2] = {Epts[i], Lzpts[j]};
                vals_even[i * N_Lz + j] = std::async(&Inversion::F_even, this, params);
                vals_odd[i * N_Lz + j] = std::async(&Inversion::F_odd, this, params);
            }
        }
        for (int i = 0; i < N_E; i++) {
            for (int j = 0; j < N_Lz; j++) {
                double val_even = vals_even[i * N_Lz + j].get();
                double val_odd = vals_odd[i * N_Lz + j].get();
                std::cout << "PSDF computed: " << Epts[i] << ", " << Lzpts[j] << std::endl;
                Fpts[(N_Lz - 1 + j) * N_E + i] = std::log(1. + val_even + val_odd);
                Fpts[(N_Lz - 1 - j) * N_E + i] = std::log(1. + val_even - val_odd);
            }
        }
    }
}

double F_even_integrand(double t, void *params) {
    inversion_params *p = (inversion_params *) params;
    Model *model = (Model *) p->model;
    InversionInterp *z2interp = (InversionInterp *) p->z2interp;
    
    std::complex<double> xi = 0.5 * p->psiEnv * (1. + std::cos(t) + 2.0i * p->h * std::sin(t));
    std::complex<double> dxi = (0.5 * p->psiEnv * (1.0i * std::sin(t) + p->h * std::cos(t)));
    std::complex<double> jac = std::pow(xi - p->E, -0.5);
    std::complex<double> R2 = std::pow(p->Lz, 2) / (2. * (xi - p->E));
    std::complex<double> z2 = model->psi_inverse(xi, p->E, p->Lz, z2interp->z2_eval(t));
    std::complex<double> drho_dpsi = model->rho_d2psi2(R2, z2, std::sqrt(R2 + z2));
    return std::real(dxi * jac * drho_dpsi);
}

double Inversion::F_even(double* params) {
    if (Inversion::verbose) std::cout << "Computing F-even: " << params[0] << ", " << params[1] << std::endl;
    
    if (params[0] == 0) return 0;
    
    double result, abserr;
    double E = params[0] * Inversion::psi0;
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

double F_odd_integrand(double t, void *params) {
    inversion_params *p = (inversion_params *) params;
    Model *model = (Model *) p->model;
    InversionInterp *z2interp = (InversionInterp *) p->z2interp;
    
    std::complex<double> xi = 0.5 * p->psiEnv * (1. + std::cos(t) + 2.0i * p->h * std::sin(t));
    std::complex<double> dxi = (0.5 * p->psiEnv * (1.0i * std::sin(t) + p->h * std::cos(t)));
    std::complex<double> jac = std::pow(xi - p->E, -0.5);
    std::complex<double> R2 = std::pow(p->Lz, 2) / (2. * (xi - p->E));
    std::complex<double> z2 = model->psi_inverse(xi, p->E, p->Lz, z2interp->z2_eval(t));
    std::complex<double> drho_dpsi = model->rho_d2psi2(R2, z2, std::sqrt(R2 + z2));
    std::complex<double> vPhi = model->v_phi(R2);
    return std::real(dxi * jac * drho_dpsi * vPhi);
}

double Inversion::F_odd(double* params) {
    if (Inversion::verbose) std::cout << "Computing F-odd: " << params[0] << ", " << params[1] << std::endl;
    
    if (!Inversion::model->is_rotating() || params[0] == 0 || params[1] == 0) return 0;
    
    double result, abserr;
    double E = params[0] * Inversion::psi0;
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

void Inversion::GSL_error_func(const char* reason, const char* file, int line, int gsl_errno) {
    std::cout << " -> GSL error #" << gsl_errno << " in line " << line << " of " << file <<": " << gsl_strerror(gsl_errno) << std::endl;
}

