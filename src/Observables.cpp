#include "Observables.hpp"

#include <fstream>
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

double minv(double a, double b) {
    return a < b ? a : b;
}

double maxv(double a, double b) {
    return a > b ? a : b;
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
    //double Rc = model->Rcirc(E * model->psi0);
    //double L = p->R * p->v * c / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    double L = p->R * p->v * c * inversion->eval_LcInv(E);
    
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
    
    return std::pow(v, 2. + p->moment) * result;
}

/**
 * Utility function for computing the numerical integral over velocity magnitude for isotropic (i.e. spherically symmetric) system
 * @param v Velocity magnitude
 * @param params Struct with the required information
 */

double rho_int_v_iso(double v, void *params) {
    struct velocity_int_params *p = (velocity_int_params *) params;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * std::pow(v, 2)) / model->psi0;
    
    return std::pow(v, 2. + p->moment) * inversion->eval_F(E, 0);
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
    
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0, 0, 0};
    
    gsl_function F;
    if (Observables::model->spherical) F.function = &rho_int_v_iso;
    else F.function = &rho_int_v;
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
    double rhoRz = Observables::rho_int(R, z, tolerance);
    
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
        p->moment = 0;
        
        double v = Observables::vMax * i / (N - 1.);
        result[2 * i] = v;
        if (v == 0 || v >= vEsc) result[2 * i + 1] = 0;
        else {
            if (Observables::model->spherical) vals[i] = std::async(rho_int_v_iso, v, p);
            else vals[i] = std::async(rho_int_v, v, p);
        }
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
    double L = 0;
    if (!model->spherical) {
        //double Rc = model->Rcirc(E * model->psi0);
        //L = p->R * vf / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
        L = p->R * vf * inversion->eval_LcInv(E);
    }
    
    return inversion->eval_F(E, L);
}

/**
 * @param v_merid Meridional velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::pv_merid_int(double v_merid, double R, double psiRz, double tolerance) {
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_merid, 0, 0};
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
    double rhoRz = Observables::rho_int(R, z, tolerance);
    
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
    double L = 0;
    if (!model->spherical) {
        //double Rc = model->Rcirc(E * model->psi0);
        //L = p->R * p->v / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
        L = p->R * p->v * inversion->eval_LcInv(E);
    }
    
    return vm * inversion->eval_F(E, L);
}

/**
 * @param v_azim Azimuthal velocity
 * @param R The value of R-coordinate
 * @param psiRz The value of gravitational potential
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::pv_azim_int(double v_azim, double R, double psiRz, double tolerance) {
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_azim, 0, 0};
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
    double rhoRz = Observables::rho_int(R, z, tolerance);
    
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
    //double Rc = model->Rcirc(E * model->psi0);
    //double L = p->R * p->vf / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    double L = p->R * p->vf * inversion->eval_LcInv(E);
    
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
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, v_rad, 0, 0};
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
    double rhoRz = Observables::rho_int(R, z, tolerance);
    
    double symmetry_factor = 1.;
    if (Observables::model->spherical) symmetry_factor = M_PI;
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        double v = Observables::vMax * (2. * i / (N - 1.) - 1.);
        result[2 * i] = v;
        if (std::abs(v) >= vEsc) result[2 * i + 1] = 0;
        else {
            if (Observables::model->spherical) vals[i] = std::async(&Observables::pv_azim_int, this, v, R, psiRz, tolerance);
            else vals[i] = std::async(&Observables::pv_rad_int, this, v, R, psiRz, tolerance);
        }
    }
    for (int i = 0; i < N; i++) {
        double v = result[2 * i];
        if (std::abs(v) < vEsc) result[2 * i + 1] = 2. * symmetry_factor / rhoRz * vals[i].get();
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

double pv_rel_wm_int(double wm, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (std::pow(wm, 2) + std::pow(p->wf, 2))) / model->psi0;
    //double Rc = model->Rcirc(E * model->psi0);
    //double L = p->R * p->wf / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    double L = p->R * p->wf * inversion->eval_LcInv(E);
    double f = inversion->eval_F(E, L);
    
    double um2 = std::pow(p->um, 2), wm2 = std::pow(wm, 2);
    double x = pow(1. - std::pow(um2 + wm2 + std::pow(p->uf - p->wf, 2) - std::pow(p->v_rel, 2), 2) / (4. * um2 * wm2), -0.5);
    return x * f;
}

double pv_rel_wf_int(double wf, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->wf = wf;
    
    double df = std::sqrt(std::pow(p->v_rel, 2) - std::pow(p->uf - wf, 2));
    double vMax = std::sqrt(2 * p->psiRz - std::pow(wf, 2));//sqrt(2. * p->P3 - pow(wf, 2) - pow(p->P6, 2) - pow(p->P7, 2));
    double vmin = maxv(df - p->um, p->um - df);
    double vmax = minv(vMax, df + p->um);
    if (vmax < vmin) return 0;

    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_wm_int;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, vmin, vmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

double pv_rel_um_int(double um, void * params) {
    
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->um = um;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * (std::pow(um, 2) + std::pow(p->uf, 2))) / model->psi0;
    //double Rc = model->Rcirc(E * model->psi0);
    //double L = p->R * p->uf / (std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc))));
    double L = p->R * p->uf * inversion->eval_LcInv(E);
    double f = inversion->eval_F(E, L);
    
    double vmin = maxv(-p->vMax, p->uf - p->v_rel);
    double vmax = minv(p->vMax, p->uf + p->v_rel);
    if (vmax < vmin) return 0;
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_wf_int;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, vmin, vmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return f * result;
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
    struct relative_velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, vMax, v_rel, 0, 0, 0, 0};
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_uf_int;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
    gsl_integration_qags(&F, -vMax, vMax, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

double pv_rel_w_int(double w, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * std::pow(w, 2)) / model->psi0;
    double f = inversion->eval_F(E, 0);
    return w * f;
}

double pv_rel_u_int(double u, void * params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->um = u;
    
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double E = (p->psiRz - 0.5 * std::pow(u, 2)) / model->psi0;
    double f = inversion->eval_F(E, 0);
    
    double vmin = std::abs(u - p->v_rel);
    double vmax = minv(p->vMax, u + p->v_rel);
    if (vmax <= vmin) return 0;
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_w_int;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, vmin, vmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return u * f * result;
}


double Observables::pv_rel_sph_int(double v_rel, double R, double psiRz, double tolerance) {
    double vMax = std::sqrt(2. * psiRz);
    struct relative_velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, vMax, v_rel, 0, 0, 0, 0};
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_u_int;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
    gsl_integration_qags(&F, 0, vMax, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return 2. * M_PI * result;
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
    double rhoRz = Observables::rho_int(R, z, tolerance);
    
    double lnvmax = std::log(2. * Observables::vMax + 1.);
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        //double v = 2. * Observables::vMax * i / (N - 1.);
        double v = std::exp(lnvmax * i / (N - 1.)) - 1.;
        result[2 * i] = v;
        if (v == 0 || v >= 2 * vEsc) result[2 * i + 1] = 0;
        else {
            if (Observables::model->spherical) vals[i] = std::async(&Observables::pv_rel_sph_int, this, v, R, psiRz, tolerance);
            else vals[i] = std::async(&Observables::pv_rel_int, this, v, R, psiRz, tolerance);
        }
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
    double rhoRz = Observables::rho_int(R, z, tolerance);
        
    struct velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, 0, 0, mom};
    gsl_function F;
    if (Observables::model->spherical) F.function = &rho_int_v_iso;
    else F.function = &rho_int_v;
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


double v_rel_mom_int(double v_rel, void *params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->v_rel = v_rel;
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_uf_int;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, -p->vMax, p->vMax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return std::pow(v_rel, 1. + p->moment) * result;
}

double v_rel_mom_sph_int(double v_rel, void *params) {
    struct relative_velocity_int_params * p = (struct relative_velocity_int_params *) params;
    p->v_rel = v_rel;
    
    double result, abserr;
    gsl_function F;
    F.function = &pv_rel_u_int;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, 0, p->vMax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return 2. * M_PI * std::pow(v_rel, 1. + p->moment) * result;
}


/**
 * @param mom Relative velocity moment (i.e. <v^{mom}>)
 * @param R The value of R-coordinate
 * @param z The value of z-coordinate
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

double Observables::v_rel_mom(int mom, double R, double z, double tolerance) {
    time_t tStart = time(NULL);
    
    double R2 = std::pow(R, 2), z2 = std::pow(z, 2);
    double psiRz = std::real(Observables::model->psi(R2, z2, std::sqrt(R2 + z2)));
    double vEsc = std::sqrt(2. * psiRz);
    double result, abserr;
    double rhoRz = Observables::rho_int(R, z, tolerance);
    
    struct relative_velocity_int_params p = {Observables::model, Observables::inversion, Observables::nIntervals, tolerance, R, psiRz, vEsc, 0, 0, 0, 0, mom};
    
    gsl_function F;
    if (Observables::model->spherical) F.function = &v_rel_mom_sph_int;
    else F.function = &v_rel_mom_int;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
    gsl_integration_qags(&F, 0, 2. * vEsc, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Relative velocity moment computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
    
    return 4. * M_PI / std::pow(rhoRz, 2) * result;
}


double occupation_int_Lz(double Lz, void *params) {
    time_t tStart = time(NULL);
    
    struct occupation_params * p = (struct occupation_params *) params;
    Inversion *inversion = (Inversion *) p->inversion;
    
    return inversion->eval_F(p->E, Lz);
}


double occupation_int_E(double E, void *params) {
    time_t tStart = time(NULL);
    
    struct occupation_params * p = (struct occupation_params *) params;
    
    Inversion *inversion = (Inversion *) p->inversion;
    Model *model = (Model *) p->model;
    p->E = E / model->psi0;
    
    double vfMax = std::sqrt(2. * (p->psiRz - E));
    
    if (model->spherical) return inversion->eval_F(p->E, 0);
    else {
        //double Rc = model->Rcirc(E);
        //double Lzc = std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc)));
        double LzcInv = inversion->eval_LcInv(p->E);
        
        double Lzmin = maxv(-p->R * vfMax * LzcInv, p->Lzmin);
        double Lzmax = minv(p->R * vfMax * LzcInv, p->Lzmax);
        if (Lzmax < Lzmin) return 0;
        
        double result, abserr;
        gsl_function F;
        F.function = &occupation_int_Lz;
        F.params = p;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
        gsl_integration_qags(&F, Lzmin, Lzmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        return result / LzcInv;
    }
}


double occupation_int_z(double z, void *params) {
    time_t tStart = time(NULL);
    
    struct occupation_params * p = (struct occupation_params *) params;
    
    Model *model = (Model *) p->model;
    
    double R2 = std::pow(p->R, 2), z2 = std::pow(p->z, 2);
    double r = std::sqrt(R2 + z2);
    double psiRz = std::real(model->psi(R2, z2, r));
    p->psiRz = psiRz;
    p->z = z;
    
    double Emin = p->Emin * model->psi0;
    double Emax = minv(p->Emax * model->psi0, psiRz);
    if (Emax < Emin) return 0;
    
    double result, abserr;
    gsl_function F;
    F.function = &occupation_int_E;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, Emin, Emax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}


double occupation_int_R(double R, void *params) {
    time_t tStart = time(NULL);
    
    struct occupation_params * p = (struct occupation_params *) params;
    p->R = R;
    
    Model *model = (Model *) p->model;
    
    double zmax = model->z_psi(p->Emin * model->psi0, R);
    double result, abserr;
    gsl_function F;
    F.function = &occupation_int_z;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, 0, zmax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

double Observables::occupation_int(double Emin, double Emax, double Lzmin, double Lzmax, double tolerance) {
    time_t tStart = time(NULL);
    
    struct occupation_params p = {Observables::model, Observables::inversion, Observables::nIntervals, Emin, Emax, Lzmin, Lzmax, tolerance, 0, 0, 0, 0};
    double Rmax = Observables::model->R_psi(Emin * Observables::model->psi0, 0);
    double result, abserr;
    gsl_function F;
    F.function = &occupation_int_R;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(Observables::nIntervals);
    gsl_integration_qags(&F, 0, Rmax, 0, tolerance, Observables::nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Occupation number (" << Emin << ", " << Emax << ", " << Lzmin << ", " << Lzmax << " -> " << result << ") computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
    */
    return 2. * std::pow(2. * M_PI, 2) * result;
}


/**
 * @param N_E Number or relative energy bins
 * @param N_Lz Number of angular momentum bins
 * @param Epts Relative energy bin edges
 * @param Lzpts Angular momentum bin edges
 * @param result An array for storing the results
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::occupation(int N_E, int N_Lz, double *Epts, double *Lzpts, double* result, double tolerance) {
    time_t tStart = time(NULL);
    
    int Npts = (N_E - 1) * (N_Lz - 1);
    /*
    double Epts[N_E + 1];
    double Lzpts[N_Lz + 1];
    
    for (int i = 0; i < N_E + 1; i++) {
        Epts[i] = 0.89 * i / (1. * N_E) + 0.1;
    }
    
    for (int i = 0; i < N_Lz + 1; i++) {
        Lzpts[i] = 2. * i / (1. * N_Lz) - 1.;
    }
    */
    
    std::vector<std::future<double>> vals(Npts);
    for (int i = 0; i < N_E - 1; i++) {
        for (int j = 0; j < N_Lz - 1; j++) {
            vals[i * (N_Lz - 1) + j] = std::async(&Observables::occupation_int, this, Epts[i], Epts[i + 1], Lzpts[j], Lzpts[j + 1], tolerance);
        }
    }
    for (int i = 0; i < Npts; i++) {
        result[i] = vals[i].get();
    }
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Occupation numbers computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
}


double dd_int_v(double v, void *params) {
    struct dd_params * p = (struct dd_params *) params;
    Model *model = (Model *) p->model;
    Inversion *inversion = (Inversion *) p->inversion;
    
    double vx = v * p->cf * p->st + p->vR_sol + p->vEarth * (0.9931 * p->cE - 0.0670 * p->sE);
    double vy = v * p->sf * p->st + p->vPhi_sol + p->vEarth * (0.1170 * p->cE + 0.4927 * p->sE);
    double vz = v * p->ct + p->vz_sol + p->vEarth * (-0.01032 * p->cE - 0.8676 * p->sE);
    
    double v2 = pow(vx, 2) + pow(vy, 2) + pow(vz, 2);
    
    double E = (p->psiR - 0.5 * v2) / model->psi0;
    if (E < 0) return 0;
    //double Rc = model->Rcirc(p->psiR - 0.5 * v2);
    //double Lzc = std::pow(Rc, 2) * std::sqrt(-2. * std::real(model->psi_dR2(std::pow(Rc, 2), 0, Rc)));
    //double Lz = p->R * vy / Lzc;
    double Lz = 0;
    if (!model->spherical) {
        Lz = p->R * vy * inversion->eval_LcInv(E);
        if (Lz < -1.) Lz = -1.;
        if (Lz > 1.) Lz = 1.;
    }
    return std::pow(v, 2. + p->power) * inversion->eval_F(E, Lz);
}


double dd_int_phi(double phi, void *params) {
    struct dd_params * p = (struct dd_params *) params;
    p->cf = std::cos(phi);
    p->sf = std::sin(phi);
    
    double cA = - (p->cf * p->st * (p->vR_sol + p->vEarth * (0.9931 * p->cE - 0.0670 * p->sE)) + p->sf * p->st * (p->vPhi_sol + p->vEarth * (0.1170 * p->cE + 0.4927 * p->sE)) + p->ct * (p->vz_sol + p->vEarth * (-0.01032 * p->cE - 0.8676 * p->sE)));
    double cB = - (pow(p->vPhi_sol, 2) + pow(p->vR_sol, 2) + pow(p->vz_sol, 2) + pow(p->vEarth, 2) + p->vEarth * (p->sE * (0.9854 * p->vPhi_sol - 0.1358 * p->vR_sol - 1.7352 * p->vz_sol) + p->cE * (0.234 * p->vPhi_sol + 1.9862 * p->vR_sol - 0.02064 * p->vz_sol - 0.00164952 * p->vEarth * p->sE)));
    double vMax = cA + sqrt(pow(cA, 2) + cB + 2 * p->psiR);
    
    if (p->vmin >= vMax) return 0;
    
    double result, abserr;
    gsl_function F;
    F.function = &dd_int_v;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, p->vmin, vMax, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return result;
}


double dd_int_theta(double theta, void *params) {
    struct dd_params * p = (struct dd_params *) params;
    p->ct = std::cos(theta);
    p->st = std::sin(theta);
    
    double result, abserr;
    gsl_function F;
    F.function = &dd_int_phi;
    F.params = p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(p->nIntervals);
    gsl_integration_qags(&F, 0, 2. * M_PI, 0, p->tolerance, p->nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    return p->st * result;
}


/**
 * @param vmin v_{min}
 * @param params DD parameters
 */

double Observables::dd_int(double vmin, dd_params params) {
    time_t tStart = time(NULL);
    
    dd_params p = {Observables::model, Observables::inversion, params.cE, params.sE, params.power, params.psiR, params.R, params.vR_sol, params.vPhi_sol, params.vz_sol, params.vEarth, vmin, 0, 0, 0, 0, params.tolerance, Observables::nIntervals};
    
    double result, abserr;
    gsl_function F;
    F.function = &dd_int_theta;
    F.params = &p;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(params.nIntervals);
    gsl_integration_qags(&F, 0, M_PI, 0, params.tolerance, params.nIntervals, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "Occupation number (" << Emin << ", " << Emax << ", " << Lzmin << ", " << Lzmax << " -> " << result << ") computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
    }
    */
    return result;
}


/**
 * @param N Number of v_{min} bins
 * @param result An array for storing the results
 * @param t Time in the year (between 0 and 1)
 * @param power Velocity power used in the itegral (set to -1 for computing g(v_{min}) and 1 for h(v_{min}))
 * @param R Radial distance
 * @param vR Peculiar radial velocity
 * @param vPhi Peculiar azimuthal velocity
 * @param vz Peculiar z velocity
 * @param vEarth Earth's circular velocity
 * @param vmax Maximum velocity to which the astrophysical factor is tabulated (can be greater then the sum of peculair motion and escape velocity!)
 * @param tolerance Relative tolerance used in performing the numerical integrals
 */

void Observables::dd(int N, double* result, double t, double power, double R, double vR, double vPhi, double vz, double vEarth, double vmax, double tolerance) {
    time_t tStart = time(NULL);
    
    double cE = std::cos(2. * M_PI * t);
    double sE = std::sin(2. * M_PI * t);
    double psiR = std::real(Observables::model->psi(std::pow(R, 2), 0, R));
    if (vPhi == -1) vPhi = std::sqrt(-2. * std::pow(R, 2) * std::real(Observables::model->psi_dR2(std::pow(R, 2), 0, R))) + 12.;
    
    dd_params params = {Observables::model, Observables::inversion, cE, sE, power, psiR, R, vR, vPhi, vz, vEarth, 0, 0, 0, 0, 0, tolerance, Observables::nIntervals};
    
    std::vector<std::future<double>> vals(N);
    for (int i = 0; i < N; i++) {
        double vmin = i * vmax / (N - 1.);
        vals[i] = std::async(&Observables::dd_int, this, vmin, params);
    }
    
    double rho_inverse = 1. / Observables::rho_int(R, 0, tolerance);
    for (int i = 0; i < N; i++) {
        result[i] = vals[i].get() * rho_inverse;
    }
    
    if (Observables::verbose) {
        double dt = difftime(time(NULL), tStart);
        std::cout << "DD astrophysical factor computed in " << (int)dt/60 << "m " << (int)dt%60 << "s!" << std::endl;
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

