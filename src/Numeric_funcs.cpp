#include "Numeric_funcs.hpp"


double psi_spheroid_0_int(double m2, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    return std::real(halo->rho(m2, 0., 0.));
}

double Numeric_funcs::psi_spheroid_0(Halo *halo, double q, double tolerance) {
    std::complex<double> e = std::sqrt(std::complex<double>(1.,0.) - std::pow(q, 2));
    std::complex<double> asin_e = std::asin(e);
    struct potential_int_params p = {halo, 0., 0., e, asin_e, q, tolerance, Numeric_funcs::nIntervals};
    
    double result, abserr;
    gsl_function F;
    F.function = &psi_spheroid_0_int;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
    gsl_integration_qagiu(&F, 0, 0, tolerance, Numeric_funcs::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    double e_factor = 1.;
    if (e != 0.) e_factor = std::real(asin_e / e);
    return 2. * M_PI * Numeric_funcs::G * q * e_factor * result;
}

double psi_spheroid_int_re(double t, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 * t + p->z2 * std::pow(std::pow(p->q, 2) - std::complex<double>(1.,0.) + 1. / t, -1);
    std::complex<double> V = p->R2 * std::pow(t, 2) + p->z2 * std::pow(std::pow(p->q, 2) - std::complex<double>(1.,0.) + 1. / t, -2);
    return std::real(halo->rho(U, 0., 0.) * V * (p->asin_e - std::asin(p->e * std::sqrt(t))) / std::pow(t, 2));
}

double psi_spheroid_int_im(double t, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 * t + p->z2 * std::pow(std::pow(p->q, 2) - std::complex<double>(1.,0.) + 1. / t, -1);
    std::complex<double> V = p->R2 * std::pow(t, 2) + p->z2 * std::pow(std::pow(p->q, 2) - std::complex<double>(1.,0.) + 1. / t, -2);
    return std::imag(halo->rho(U, 0., 0.) * V * (p->asin_e - std::asin(p->e * std::sqrt(t))) / std::pow(t, 2));
}

double psi_spheroid_int_q1_re(double t, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = (p->R2 + p->z2) * t;
    std::complex<double> V = U * t;
    return std::real(halo->rho(U, 0., 0.) * V * (1. - std::sqrt(t)) / std::pow(t, 2));
}

double psi_spheroid_int_q1_im(double t, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = (p->R2 + p->z2) * t;
    std::complex<double> V = U * t;
    return std::imag(halo->rho(U, 0., 0.) * V * (1. - std::sqrt(t)) / std::pow(t, 2));
}

std::complex<double> Numeric_funcs::psi_spheroid(Halo* halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance) {
    std::complex<double> e = 1.;
    if (q != 1.) e = std::sqrt(std::complex<double>(1.,0.) - std::pow(q, 2));
    std::complex<double> asin_e = std::asin(e);
    struct potential_int_params p = {halo, R2, z2, e, asin_e, q, tolerance, Numeric_funcs::nIntervals};
    
    double result_re, result_im, abserr_re, abserr_im;
    gsl_function F;
    if (q == 1.) F.function = &psi_spheroid_int_q1_re;
    else F.function = &psi_spheroid_int_re;
    F.params = &p;
    gsl_integration_workspace *workspace_re = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
    gsl_integration_qag(&F, 0., 1., 0., tolerance, Numeric_funcs::nIntervals, 6, workspace_re, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace_re);
    
    if (std::imag(R2) == 0 && std::imag(z2) == 0) result_im = 0;
    else {
        if (q == 1.) F.function = &psi_spheroid_int_q1_im;
        else F.function = &psi_spheroid_int_im;
        F.params = &p;
        gsl_integration_workspace *workspace_im = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
        gsl_integration_qag(&F, 0., 1., 0., tolerance, Numeric_funcs::nIntervals, 6, workspace_im, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace_im);
    }
    return 2. * M_PI * Numeric_funcs::G * q / e * std::complex<double>(result_re, result_im);
}


double psi_spheroid_dR2_int_re(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::real(halo->rho(U, 0., 0.) / (std::pow(1. + u, 2) * std::sqrt(std::pow(p->q, 2) + u)));
}

double psi_spheroid_dR2_int_im(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::imag(halo->rho(U, 0., 0.) / (std::pow(1. + u, 2) * std::sqrt(std::pow(p->q, 2) + u)));
}

std::complex<double> Numeric_funcs::psi_spheroid_dR2(Halo* halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance) {
    struct potential_int_params p = {halo, R2, z2, 0., 0., q, tolerance, Numeric_funcs::nIntervals};
    
    double result_re, result_im, abserr_re, abserr_im;
    gsl_function F;
    F.function = &psi_spheroid_dR2_int_re;
    F.params = &p;
    gsl_integration_workspace *workspace_re = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
    gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_re, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace_re);
    
    if (std::imag(R2) == 0 && std::imag(z2) == 0) result_im = 0;
    else {
        F.function = &psi_spheroid_dR2_int_im;
        F.params = &p;
        gsl_integration_workspace *workspace_im = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
        gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_im, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace_im);
    }
    return - M_PI * Numeric_funcs::G * q * std::complex<double>(result_re, result_im);
}

double psi_spheroid_dz2_int_re(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::real(halo->rho(U, 0., 0.) / ((1. + u) * std::pow(std::pow(p->q, 2) + u, 1.5)));
}

double psi_spheroid_dz2_int_im(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::imag(halo->rho(U, 0., 0.) / ((1. + u) * std::pow(std::pow(p->q, 2) + u, 1.5)));
}

std::complex<double> Numeric_funcs::psi_spheroid_dz2(Halo* halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance) {
    struct potential_int_params p = {halo, R2, z2, 0., 0., q, tolerance, Numeric_funcs::nIntervals};
    
    double result_re, result_im, abserr_re, abserr_im;
    gsl_function F;
    F.function = &psi_spheroid_dz2_int_re;
    F.params = &p;
    gsl_integration_workspace *workspace_re = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
    gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_re, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace_re);
    
    if (std::imag(R2) == 0 && std::imag(z2) == 0) result_im = 0;
    else {
        F.function = &psi_spheroid_dz2_int_im;
        F.params = &p;
        gsl_integration_workspace *workspace_im = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
        gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_im, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace_im);
    }
    return - M_PI * Numeric_funcs::G * q * std::complex<double>(result_re, result_im);
}

double psi_spheroid_d2z2_int_re(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::real(halo->rho_dz2(U, 0., 0.) / ((1. + u) * std::pow(std::pow(p->q, 2) + u, 2.5)));
}

double psi_spheroid_d2z2_int_im(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::imag(halo->rho_dz2(U, 0., 0.) / ((1. + u) * std::pow(std::pow(p->q, 2) + u, 2.5)));
}

std::complex<double> Numeric_funcs::psi_spheroid_d2z2(Halo* halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance) {
    struct potential_int_params p = {halo, R2, z2, 0., 0., q, tolerance, Numeric_funcs::nIntervals};
    
    double result_re, result_im, abserr_re, abserr_im;
    gsl_function F;
    F.function = &psi_spheroid_d2z2_int_re;
    F.params = &p;
    gsl_integration_workspace *workspace_re = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
    gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_re, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace_re);
    
    if (std::imag(R2) == 0 && std::imag(z2) == 0) result_im = 0;
    else {
        F.function = &psi_spheroid_d2z2_int_im;
        F.params = &p;
        gsl_integration_workspace *workspace_im = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
        gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_im, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace_im);
    }
    return - M_PI * Numeric_funcs::G * std::pow(q, 3) * std::complex<double>(result_re, result_im);
}

double psi_spheroid_d2R2z2_int_re(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::real(halo->rho_dz2(U, 0., 0.) / (std::pow(1. + u, 2) * std::pow(std::pow(p->q, 2) + u, 1.5)));
}

double psi_spheroid_d2R2z2_int_im(double u, void *params) {
    struct potential_int_params *p = (potential_int_params *) params;
    Halo *halo = (Halo *) p->halo;
    
    std::complex<double> U = p->R2 / (1. + u) + p->z2 / (std::pow(p->q, 2) + u);
    return std::imag(halo->rho_dz2(U, 0., 0.) / (std::pow(1. + u, 2) * std::pow(std::pow(p->q, 2) + u, 1.5)));
}

std::complex<double> Numeric_funcs::psi_spheroid_d2R2z2(Halo* halo, std::complex<double> R2, std::complex<double> z2, double q, double tolerance) {
    struct potential_int_params p = {halo, R2, z2, 0., 0., q, tolerance, Numeric_funcs::nIntervals};
    
    double result_re, result_im, abserr_re, abserr_im;
    gsl_function F;
    F.function = &psi_spheroid_d2z2_int_re;
    F.params = &p;
    gsl_integration_workspace *workspace_re = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
    gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_re, &result_re, &abserr_re);
    gsl_integration_workspace_free(workspace_re);
    
    if (std::imag(R2) == 0 && std::imag(z2) == 0) result_im = 0;
    else {
        F.function = &psi_spheroid_d2R2z2_int_im;
        F.params = &p;
        gsl_integration_workspace *workspace_im = gsl_integration_workspace_alloc(Numeric_funcs::nIntervals);
        gsl_integration_qagiu(&F, 0., 0., tolerance, Numeric_funcs::nIntervals, workspace_im, &result_im, &abserr_im);
        gsl_integration_workspace_free(workspace_im);
    }
    return - M_PI * Numeric_funcs::G * std::pow(q, 3) * std::complex<double>(result_re, result_im);
}


std::complex<double> Numeric_funcs::beta_cf(std::complex<double> x, double a, double b, int N, int iterations) {
    if (N > iterations) return 0;
    std::complex<double> dn;
    int m = N / 2;
    if (N % 2 == 0) dn = m * (b - m) * x / ((a + 2. * m - 1.) * (a + 2. * m));
    else dn = - (a + m) * (a + b + m) * x / ((a + 2. * m) * (a + 2. * m + 1.));
    return dn / (1. + Numeric_funcs::beta_cf(x, a, b, N + 1, iterations));
}

std::complex<double> Numeric_funcs::psi_gNFW(halo_6p *params, std::complex<double> r, int iterations) {
    if (r == 0.) return 4. * M_PI * Numeric_funcs::G * params->rho_s * std::pow(params->r_s, 2) / (params->gamma - 2);
    double a = 3. - params->gamma, b = params->gamma - 1;
    std::complex<double> x = r / params->r_s;
    std::complex<double> beta_val = std::pow(-x, a) * std::pow(1. + x, b) / a / (1. + beta_cf(-x, a, b, 1, iterations));
    return -4. * M_PI * Numeric_funcs::G * params->rho_s * std::pow(params->r_s, 2) * (x + std::pow(std::complex<double>(-1.,0.), params->gamma) * beta_val) / (x * (params->gamma - 2));
}

std::complex<double> Numeric_funcs::psi_gNFW_dr2(halo_6p *params, std::complex<double> r, int iterations) {
    if (r == 0.) return 0.;
    double a = 3. - params->gamma, b = params->gamma - 2;
    std::complex<double> x = r / params->r_s;
    std::complex<double> beta_val = std::pow(-x, a) * std::pow(1. + x, b) / a / (1. + beta_cf(-x, a, b, 1, iterations));
    return -2. * M_PI * Numeric_funcs::G * params->rho_s * std::pow(params->r_s, 3) * std::pow(std::complex<double>(-1.,0.), params->gamma + 1.) * beta_val * std::pow(r, -3);
}

std::complex<double> Numeric_funcs::psi_gNFW_d2r2(Halo *halo, std::complex<double> r) {
    if (r == 0.) return 0.;
    std::complex<double> r2 = std::pow(r, 2);
    return -(4. * M_PI * Numeric_funcs::G * halo->rho(r2, 0, r) + 6. * halo->psi_dR2(r2, 0, r)) / (4. * r2);
}

/**
 * @param reason Reason for raising the error
 * @param file Name of the file in which the error occured
 * @param line Line in which the error occured
 * @param gsl_errno GSL error code
 */

void Numeric_funcs::GSL_error_func(const char* reason, const char* file, int line, int gsl_errno) {
    std::cout << " -> GSL error #" << gsl_errno << " in line " << line << " of " << file <<": " << gsl_strerror(gsl_errno) << std::endl;
}

