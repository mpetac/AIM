#include "Observables.hpp"

/**
 * @param model Galactic model
 * @param inv Instance of the PSDF inversion class
 */

Observables::Observables(Model *model, Inversion *inversion) {
    gsl_set_error_handler(&Observables::GSL_error_func);
    Observables::model = model;
    Observables::inversion = inversion;
}



/**
 * Utility function for computing the numerical integral over velocity angle
 * @param c Cosine of the angle
 * @param params Struct with the required information
 */

double rho_int_c(double c, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * std::pow(p->v,2)) / model->psi0;
    double Rc = model->Rcirc(E * model->psi0);
    double L = p->R * p->v * c / (std::pow(Rc, 2) * sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    
    return inversion->eval_F(E, L);
}

/**
 * Utility function for computing the numerical integral over velocity magnitude
 * @param v Velocity magnitude
 * @param params Struct with the required information
 */

double rho_int_v(double v, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    p->v = v;
    
    double result, abserr; 
    
    gsl_function F;
    F.function = &rho_int_c;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, 0, 1, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return std::pow(v, 2) * result;
}

/**
 * @param R Value of the R-coordinate
 * @param z Value of the z-coordinate
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::rho_int(double R, double z, double tolerance) {
    double result, abserr;
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0};
    
    gsl_function F;
    F.function = &rho_int_v;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, std::sqrt(2. * psiRz), 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return 4. * M_PI * result;
}

/**
 * @param N Number of points in which DM density will be computed
 * @param Rpts The values of R-coordinates in which DM density will be computed
 * @param zpts The values of z-coordinates in which DM density will be computed
 * @param result An array for stroing the DM density
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::rho(int N, double* Rpts, double* zpts, double* result, double tolerance) {
    for (int i = 0; i < N; i++) {
        result[i] = Observables::rho_int(Rpts[i], zpts[i], tolerance);
    }
}

/**
 * @param N Number of points in which the velocity distribution will be computed
 * @param R The value of R-coordinate
 * @param z The value of z-coordinate
 * @param result An array for storing the velocity probability distribution
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::pv_mag(int N, double R, double z, double* result, double tolerance) {
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    for (int i = 0; i < N; i++) {
        struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0};
        double v = vEsc * i / (N - 1.);
        result[2 * i] = v;
        if (v == 0 || v == vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 4 * M_PI / rhoRz * rho_int_v(v, &p);
    }
}

/**
 * Utility function for computing the numerical integral over azimuthal velocity
 * @param vf Azimuthal velocity
 * @param params Struct with the required information
 */

double pv_merid_int_vf(double vf, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (std::pow(p->v, 2) + std::pow(vf, 2))) / model->psi0;
    double Rc = model->Rcirc(E * model->psi0);
    double L = p->R * vf / (std::pow(Rc, 2) * sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    
    return inversion->eval_F(E, L);
}

/**
 * @param v_merid Meridional velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::pv_merid_int(double v_merid, double R, double psiRz, double tolerance) {
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_merid};
    double vMax = std::sqrt(2. * psiRz - std::pow(v_merid, 2));
    double result, abserr;
    gsl_function F;
    F.function = &pv_merid_int_vf;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, -vMax, vMax, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

/**
 * @param N Number of points in which the velocity distribution will be computed
 * @param R The value of R-coordinate
 * @param z The value of z-coordinate
 * @param result An array for storing the velocity probability distribution
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::pv_merid(int N, double R, double z, double* result, double tolerance) {
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    for (int i = 0; i < N; i++) {
        double v_merid = vEsc * i / (N - 1.);
        result[2 * i] = v_merid;
        if (v_merid == 0 || v_merid == vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 2 * M_PI * v_merid / rhoRz * Observables::pv_merid_int(v_merid, R, psiRz, tolerance);
    }
}

/**
 * Utility function for computing the numerical integral over meridional velocity
 * @param vf Azimuthal velocity
 * @param params Struct with the required information
 */

double pv_azim_int_vm(double vm, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (std::pow(p->v, 2) + std::pow(vm, 2))) / model->psi0;
    double Rc = model->Rcirc(E * model->psi0);
    double L = p->R * p->v / (std::pow(Rc, 2) * sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    
    return vm * inversion->eval_F(E, L);
}

/**
 * @param v_azim Azimuthal velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::pv_azim_int(double v_azim, double R, double psiRz, double tolerance) {
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_azim};
    double vMax = std::sqrt(2. * psiRz - std::pow(v_azim, 2));
    double result, abserr;
    gsl_function F;
    F.function = &pv_azim_int_vm;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nIntervals);
    gsl_integration_qags(&F, 0, vMax, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

/**
 * @param N Number of points in which the velocity distribution will be computed
 * @param R The value of R-coordinate
 * @param z The value of z-coordinate
 * @param result An array for storing the velocity probability distribution
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::pv_azim(int N, double R, double z, double* result, double tolerance) {
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    for (int i = 0; i < N; i++) {
        double v_azim = vEsc * (2. * i / (N - 1.) - 1.);
        result[2 * i] = v_azim;
        if (std::abs(v_azim) >= vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 2. * M_PI / rhoRz * Observables::pv_azim_int(v_azim, R, psiRz, tolerance);
    }
}


/**
 * @param reason Reason for raising the error
 * @param file Name of the file in which the error occured
 * @param line Line in which the error occured
 * @param gsl_errno GSL error code
 */

void Observables::GSL_error_func(const char* reason, const char* file, int line, int gsl_errno) {
    std::cout << " -> GSL error #" << gsl_errno << " in line " << line << " of " << file <<": " << gsl_strerror(gsl_errno) << std::endl;
}

