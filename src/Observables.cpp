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
 * @param N Number of points in which DM density will be computed
 * @param Rpts The values of R-coordinates in which DM density will be computed
 * @param Rpts The values of z-coordinates in which DM density will be computed
 * @param result An array in which the resulting DM densities will be stored
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::rho(int N, double* Rpts, double* zpts, double* result, double tolerance) {
    for (int i = 0; i < N; i++) {
        result[i] = rho_int(Rpts[i], zpts[i], tolerance);
    }
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
    
    double E = (p->psiRz - 0.5 * pow(p->v, 2)) / model->psi0;
    
    double Rc = model->Rcirc(E * model->psi0);
    double L = p->R * p->v * c / (std::pow(Rc, 2) * sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    //double L = p->R * p->v * c * inversion->eval_LcI(E);
    
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
 * @param reason Reason for raising the error
 * @param file Name of the file in which the error occured
 * @param line Line in which the error occured
 * @param gsl_errno GSL error code
 */

void Observables::GSL_error_func(const char* reason, const char* file, int line, int gsl_errno) {
    std::cout << " -> GSL error #" << gsl_errno << " in line " << line << " of " << file <<": " << gsl_strerror(gsl_errno) << std::endl;
    std::cout << " -> " << reason << std::endl;
}

