#include "Observables.hpp"

/**
 * @param model Galactic model
 * @param inv Instance of the PSDF inversion class
 */

Observables::Observables(Model *model, Inversion *inversion, bool verbose) {
    Observables::model = model;
    Observables::inversion = inversion;
    Observables::verbose = verbose;
    Observables::vMax = std::sqrt(2. * std::real(model->psi0));
    
    if (verbose) gsl_set_error_handler(&Observables::GSL_error_func);
    else gsl_set_error_handler(&Observables::GSL_error_func_silent);
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
    
    return std::pow(v, 2 + p->vf) * result;
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
    
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0, 0};
    
    gsl_function F;
    F.function = &rho_int_v;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
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
    time_t tStart = time(NULL);
    
    //for (int i = 0; i < N; i++)  result[i] = Observables::rho_int(Rpts[i], zpts[i], tolerance);
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) vals[i] = std::async(&Observables::rho_int, this, Rpts[i], zpts[i], tolerance);
    for (int i = 0; i < N; i++) result[i] = vals[i].get();
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Numerical density profile computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
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
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    //double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double rhoRz = Observables::rho_int(R2, z2, tolerance);
    
    
    /*
    for (int i = 0; i < N; i++) {
        struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0, 0};
        double v = Observables::vMax * i / (N - 1.);
        result[2 * i] = v;
        if (v == 0 || v >= vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 4. * M_PI / rhoRz * rho_int_v(v, &p);
    }
    */
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        struct velocity_int_params *p = new velocity_int_params();
        p->model = Observables::model;
        p->inversion = Observables::inversion;
        p->nIntervals = Observables::nIntervals;
        p->tolerance = tolerance;
        p->R = R;
        p->psiRz = psiRz;
        p->v = 0;
        p->vf = 0;
        
        double v = Observables::vMax * i / (N - 1.);
        result[2 * i] = v;
        if (v == 0 || v >= vEsc) result[2 * i + 1] = 0;
        else vals[i] = std::async(rho_int_v, v, p);
    }
    for (int i = 0; i < N; i++) {
        double v = result[2 * i];
        if (v > 0 && v < vEsc) result[2 * i + 1] = 4. * M_PI / rhoRz * vals[i].get();
    }
    
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Velocity magnitude distribution computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
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
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_merid, 0};
    double vMax = std::sqrt(2. * psiRz - std::pow(v_merid, 2));
    double result, abserr;
    gsl_function F;
    F.function = &pv_merid_int_vf;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
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
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    //double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double rhoRz = Observables::rho_int(R2, z2, tolerance);
    
    /*
    for (int i = 0; i < N; i++) {
        double v_merid = Observables::vMax * i / (N - 1.);
        result[2 * i] = v_merid;
        if (v_merid == 0 || v_merid >= vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 2 * M_PI * v_merid / rhoRz * Observables::pv_merid_int(v_merid, R, psiRz, tolerance);
    }
    */
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        double v = Observables::vMax * i / (N - 1.);
        result[2 * i] = v;
        if (v == 0 || v >= vEsc) result[2 * i + 1] = 0;
        else vals[i] = std::async(&Observables::pv_merid_int, this, v, R, psiRz, tolerance);
    }
    for (int i = 0; i < N; i++) {
        double v = result[2 * i];
        if (v > 0 && v < vEsc) result[2 * i + 1] = 2. * M_PI * v / rhoRz * vals[i].get();
    }
    
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Meridional velocity distribution computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
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
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_azim, 0};
    double vMax = std::sqrt(2. * psiRz - std::pow(v_azim, 2));
    double result, abserr;
    gsl_function F;
    F.function = &pv_azim_int_vm;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
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
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    //double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double rhoRz = Observables::rho_int(R2, z2, tolerance);
    
    /*
    for (int i = 0; i < N; i++) {
        double v_azim = Observables::vMax * (2. * i / (N - 1.) - 1.);
        result[2 * i] = v_azim;
        if (std::abs(v_azim) >= vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 2. * M_PI / rhoRz * Observables::pv_azim_int(v_azim, R, psiRz, tolerance);
    }
    */
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        double v = Observables::vMax * (2. * i / (N - 1.) - 1.);
        result[2 * i] = v;
        if (std::abs(v) >= vEsc) result[2 * i + 1] = 0;
        else vals[i] = std::async(&Observables::pv_azim_int, this, v, R, psiRz, tolerance);
    }
    for (int i = 0; i < N; i++) {
        double v = result[2 * i];
        if (std::abs(v) < vEsc) result[2 * i + 1] = 2. * M_PI / rhoRz * vals[i].get();
    }
    
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Azimuthal velocity distribution computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
}


/**
 * Utility function for computing the numerical integral over meridional velocity
 * @param vf Azimuthal velocity
 * @param params Struct with the required information
 */

double pv_rad_int_vz(double vz, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (std::pow(p->v, 2) + std::pow(p->vf, 2) + std::pow(vz, 2))) / model->psi0;
    double Rc = model->Rcirc(E * model->psi0);
    double L = p->R * p->vf / (std::pow(Rc, 2) * sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    
    return inversion->eval_F(E, L);
}

/**
 * @param v_azim Azimuthal velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double pv_rad_int_vf(double vf, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    p->vf = vf;
    
    double vMax = std::sqrt(2. * p->psiRz - std::pow(p->v, 2) - std::pow(vf, 2));
    double result, abserr;
    gsl_function F;
    F.function = &pv_rad_int_vz;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, 0, vMax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}


/**
 * @param v_azim Azimuthal velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::pv_rad_int(double v_rad, double R, double psiRz, double tolerance) {
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_rad, 0};
    double vMax = std::sqrt(2. * psiRz - std::pow(v_rad, 2));
    double result, abserr;
    gsl_function F;
    F.function = &pv_rad_int_vf;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
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

void Observables::pv_rad(int N, double R, double z, double* result, double tolerance) {
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    //double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double rhoRz = Observables::rho_int(R2, z2, tolerance);
    
    /*
    for (int i = 0; i < N; i++) {
        double v_rad = Observables::vMax * (2. * i / (N - 1.) - 1.);
        result[2 * i] = v_rad;
        if (std::abs(v_rad) >= vEsc) result[2 * i + 1] = 0;
        else result[2 * i + 1] = 2. / rhoRz * Observables::pv_rad_int(v_rad, R, psiRz, tolerance);
    }
    */
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        double v = Observables::vMax * (2. * i / (N - 1.) - 1.);
        result[2 * i] = v;
        if (std::abs(v) >= vEsc) result[2 * i + 1] = 0;
        else vals[i] = std::async(&Observables::pv_rad_int, this, v, R, psiRz, tolerance);
    }
    for (int i = 0; i < N; i++) {
        double v = result[2 * i];
        if (std::abs(v) < vEsc) result[2 * i + 1] = 2. / rhoRz * vals[i].get();
    }
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Radial velocity distribution computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
}


/**
 * @param v_azim Azimuthal velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */


double minv(double a, double b) {
    return a < b ? a : b;
}

double maxv(double a, double b) {
    return a > b ? a : b;
}

double pv_rel_wm_int(double wm, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (pow(wm, 2) + pow(p->wf, 2))) / model->psi0;
    //double Rc = model->Rcirc(E * model->psi0);
    //double L = p->R * p->wf / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(Rc * Rc, 0, Rc))));
    double L = maxv(minv(p->R * p->wf * inversion->eval_LcInv(E), 1), -1);
    double f = inversion->eval_F(E, L);
    
    double um2 = pow(p->um, 2), wm2 = pow(wm, 2);
    double x = pow(1. - pow(um2 + wm2 + pow(p->uf - p->wf, 2) - pow(p->v_rel, 2), 2) / (4. * um2 * wm2), -0.5);
    return x * f;
}

double pv_rel_wf_int(double wf, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->wf = wf;
    
    double df = std::sqrt(std::pow(p->v_rel, 2) - std::pow(p->uf - wf, 2));
    double vMax = std::sqrt(2 * p->psiRz - pow(wf, 2));//sqrt(2. * p->P3 - pow(wf, 2) - pow(p->P6, 2) - pow(p->P7, 2));
    double vmin = maxv(df - p->um, p->um - df);
    double vmax = minv(vMax, df + p->um);
    
    if (vmin < vmax) {
        double result, abserr;
        gsl_function F;
        F.function = &pv_rel_wm_int;
        F.params = p;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
        gsl_integration_qags(&F, vmin, vmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return result;
    } else {
        //printf("vmin > vmax in wf_integrand (vmin: %g, vmax: %g, vesc: %g)\n", vmin, vmax, vMax);
        return 0;
    }
}

double pv_rel_um_int(double um, void * params) {
    
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->um = um;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (pow(um, 2) + pow(p->uf, 2))) / model->psi0;
    //double Rc = model->Rcirc(E * model->psi0);
    //double L = p->R * p->uf / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(Rc * Rc, 0, Rc))));
    double L = maxv(minv(p->R * p->wf * inversion->eval_LcInv(E), 1), -1);
    double f = inversion->eval_F(E, L);
    
    
    double vmin = maxv(-p->vMax, p->uf - p->v_rel);
    double vmax = minv(p->vMax, p->uf + p->v_rel);
    if (vmin < vmax) {
        double result, abserr;
        gsl_function F;
        F.function = &pv_rel_wf_int;
        F.params = p;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
        gsl_integration_qags(&F, vmin, vmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return f * result;
    } else {
        //printf("vmin > vmax in um_integrand\n");
        return 0;
    }
}


double pv_rel_uf_int(double uf, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->uf = uf;
    double vMax = std::sqrt(2. * p->psiRz - std::pow(uf, 2));
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_um_int;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, 0, vMax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}


double Observables::pv_rel_int(double v_rel, double R, double psiRz, double tolerance) {
    double vMax = std::sqrt(2. * psiRz);
    struct relative_velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, vMax, v_rel, 0, 0, 0};
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_uf_int;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
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

void Observables::pv_rel(int N, double R, double z, double* result, double tolerance) {
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    //double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double rhoRz = Observables::rho_int(R2, z2, tolerance);
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        double v = 2. * Observables::vMax * i / (N - 1.);
        result[2 * i] = v;
        if (v == 0 || v >= 2 * vEsc) result[2 * i + 1] = 0;
        else vals[i] = std::async(&Observables::pv_rel_int, this, v, R, psiRz, tolerance);
    }
    for (int i = 0; i < N; i++) {
        double v = result[2 * i];
        if (v > 0 && v < 2 * vEsc) result[2 * i + 1] = 4. * M_PI * v / std::pow(rhoRz, 2) * vals[i].get();
    }
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Relative velocity distribution computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
}


/**
 * @param mom Velocity moment (i.e. <v^{mom}>)
 * @param R The value of R-coordinate
 * @param z The value of z-coordinate
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::v_mom(int mom, double R, double z, double tolerance) {
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    double result, abserr;
    //double rhoRz = std::real(Observables::model->rho(R2, z2, std::sqrt(R2 + z2)));
    double rhoRz = Observables::rho_int(R2, z2, tolerance);
        
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0, 1. * mom};
    //    rho_int_v(v, &p);
    gsl_function F;
    F.function = &rho_int_v;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
    gsl_integration_qags(&F, 0, vEsc, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Velocity moment computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
    return 4. * M_PI / rhoRz * result; 
}


double Observables::occupation_int(double Emin, double Emax, double Lzmin, double Lzmax, double tolerance) {
    
}


void Observables::occupation(int N_E, int N_Lz, double* result, double tolerance) {
    double Epts[N_E];
    double Lzpts[N_Lz];
    
    for (int i = 0; i < N_E; i++) {
        Epts[i] = 1. * i / (N_E - 1.);
    }
    
    for (int i = 0; i < N_Lz; i++) {
        Lzpts[i] = 2. * i / (N_Lz - 1.) - 1.;
    }
    
    for (int i = 0; i < N_E - 1; i++) {
        for (int j = 0; j < N_Lz - 1; j++) {
            result[i * N_Lz + j] = Observables::occupation_int(Epts[i], Epts[i + 1], Lzpts[j], Lzpts[j + 1], tolerance);
        }
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

/**
 * @param reason Reason for raising the error
 * @param file Name of the file in which the error occured
 * @param line Line in which the error occured
 * @param gsl_errno GSL error code
 */

void Observables::GSL_error_func_silent(const char* reason, const char* file, int line, int gsl_errno) {
}

