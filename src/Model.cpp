#include "Model.hpp"

Model::Model(Halo* halo_model, Baryons* baryons_model, bool v) {
    Model::halo = halo_model;
    Model::baryons = baryons_model;
    Model::verbose = v;
    Model::psi0 = std::real(Model::psi(0, 0, 0));
    Model::spherical = Model::psi(1., 0, 1.) == Model::psi(0, 1., 1.);
    Model::rotating = Model::halo->is_rotating();
}


/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi(R2, z2, r) + Model::baryons->psi(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_dR2(R2, z2, r) + Model::baryons->psi_dR2(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_dz2(R2, z2, r) + Model::baryons->psi_dz2(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_d2R2z2(R2, z2, r) + Model::baryons->psi_d2R2z2(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_d2z2(R2, z2, r) + Model::baryons->psi_d2z2(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->rho(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->rho_dz2(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->rho_d2z2(R2, z2, r);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::rho_d2psi2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> dpsi_dz2 = Model::psi_dz2(R2, z2, r);
    return Model::rho_d2z2(R2, z2, r) * std::pow(dpsi_dz2, -2) - Model::rho_dz2(R2, z2, r) * Model::psi_d2z2(R2, z2, r) * std::pow(dpsi_dz2, -3);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 * @param r Value of radial coordinate
 */

std::complex<double> Model::rho_vphi_d2psi2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> dpsi_dz2 = Model::psi_dz2(R2, z2, r);
    std::complex<double> rho = Model::rho(R2, z2, r);
    std::complex<double> drho_dz2 = Model::rho_dz2(R2, z2, r);
    std::complex<double> vphi = Model::v_phi(R2, z2);
    std::complex<double> dvphi_dz2 = Model::v_phi_dz2(R2, z2);
    return (Model::rho_d2z2(R2, z2, r) * vphi + 2. * drho_dz2 * dvphi_dz2 + rho * Model::v_phi_d2z2(R2, z2)) * std::pow(dpsi_dz2, -2) - (drho_dz2 * vphi + rho * dvphi_dz2) * Model::psi_d2z2(R2, z2, r) * std::pow(dpsi_dz2, -3);
}


/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 */

std::complex< double > Model::v_phi(std::complex< double > R2, std::complex< double > z2) {
    return Model::halo->v_phi(R2, z2);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 */

std::complex<double> Model::v_phi_dz2(std::complex<double> R2, std::complex<double> z2) {
    return Model::halo->v_phi_dz2(R2, z2);
}

/** 
 * @param R2 Value of the R-coordinate squared
 * @param z2 Value of the z-coordinate squared
 */

std::complex<double> Model::v_phi_d2z2(std::complex<double> R2, std::complex<double> z2) {
    return Model::halo->v_phi_d2z2(R2, z2);
}

/** 
 * Computes the relative difference between psi(R2,z2) and xi
 * @param v Vector of complex and imaginary value of z2 that is being searched for
 * @param p Struct containing the remianing relevant information
 */

double psi_inverse_diff(const gsl_vector* v, void* p) {
    psi_inverse_params *params = (psi_inverse_params *) p;
    Model *model = (Model *) params->model;
    std::complex<double> z2(gsl_vector_get(v, 0), gsl_vector_get(v, 1));
    std::complex<double> r = std::sqrt(params->R2 + z2);
    std::complex<double> psi_val = model->psi(params->R2, z2, r);
    //printf("Psi(%g+%gI, %g+%gI): %g+%gI -> diff: %.32g\n", std::real(params->R2), std::imag(params->R2), std::real(z2), std::imag(z2), std::real(psi_val), std::imag(psi_val), std::abs((params->xi - psi_val) / params->xi));
    return std::abs((params->xi - psi_val) / psi_val);
}

/** 
 * @param xi Target value of total gravitational potential
 * @param E Relative energy
 * @param Lz Angular momentum around z-axis
 * @param z0 Inital guess for z2
 * @param tolerance Requested relative tolerance
 * @param limit Maximum number of iterations
 */

std::complex<double> Model::psi_inverse(std::complex<double> xi, double E, double Lz, std::complex<double> z0, double tolerance, int limit) {
    int i = 0;
    double diff = 1;
    //if (std::abs(z0) > 1e4) limit = 1000;
    
    psi_inverse_params params = {this, xi, std::pow(Lz, 2) / (2. * (xi - E))};
    
    gsl_multimin_function minf;
    minf.n = 2;
    minf.params = &params;
    minf.f = psi_inverse_diff;
    
    gsl_vector *x, *step;
    x = gsl_vector_alloc(2);
    step = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, std::real(z0));
    gsl_vector_set(x, 1, std::imag(z0));
    gsl_vector_set(step, 0, std::fabs(std::real(z0) + 1e-1) / 10.);
    gsl_vector_set(step, 1, std::fabs(std::imag(z0) + 1e-1) / 10.);
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 2);
    gsl_multimin_fminimizer_set(s, &minf, x, step);
    
    while (diff > tolerance && i < limit) {
        i++;
        gsl_multimin_fminimizer_iterate(s);
        diff = s->fval;
    }
    
    std::complex<double> result(gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
    gsl_vector_free(x);
    gsl_vector_free(step);
    gsl_multimin_fminimizer_free(s);
    
    /**/
    if (diff > tolerance) {
        //std::cout << "Inversion failed! xi: " << xi << ", E: " << E << ", Lz: " << Lz << ", iter: " << i << ", delta: " << s->fval << ", z0: " << z0 << " -> " << result << std::endl;
        //printf("Inversion failed! (xi: %.15g+%.15gI, E: %.15g, L: %.15g, iter: %d, delta: %g, z0: %g+%gI -> %g+%gI)\n", std::real(xi), std::imag(xi), E, Lz, i, s->fval, std::real(z0), std::imag(z0), std::real(result), std::imag(result));
    }
    
    return result;
}

/** 
 * @param E Relative energy
 * @param tolerance Requested relative tolerance
 * @param limit Maximum number of iterations
 */

double Model::Rcirc(double E, double tolerance, int limit) {
    if (E == 0) return 1e306;
    int i = 0;
    double Rc = 0, dR = 1.;
    while(i < limit) {
        i++;
        double R = Rc + dR;
        double R2 = std::pow(R, 2);
        double E_R = std::real(Model::psi(R2, 0, R) + R2 * Model::psi_dR2(R2, 0, R));
        
        //printf("%3.d E = %g: R = %g, E_R = %g (Rc = %g, dR = %g)\n", i, E, R, E_R, Rc, dR);
        
        if (E > E_R) dR /= 2.;
        else {
            Rc = R;
            dR *= 2.;
        }
        
        if (std::abs(E_R - E) < tolerance * E) {
            Rc = R;
            break;
        }
    }
    
    if (Model::verbose && i >= limit) std::cout << "Rc inversion did not converge! E = " << E << "-> Rc = " << Rc << ", E(Rc) = " << std::real(Model::psi(Rc * Rc, 0, Rc) + (Rc * Rc) * Model::psi_dR2(Rc * Rc, 0, Rc)) << std::endl;
    return Rc;
}


double Model::R_psi(double psi, double z, double tolerance, int limit) {
    if (psi == 0) return 1e306;
    
    double z2 = std::pow(z, 2);
    double R0 = 0, dR = 1.;
    int i = 0;
    while(i < limit) {
        i++;
        double R = R0 + dR;
        double R2 = std::pow(R, 2);
        double r = std::sqrt(R2 + z2);
        double psi_R = std::real(Model::psi(R2, z2, r));
        
        if (psi > psi_R) dR /= 2.;
        else {
            R0 = R;
            dR *= 2.;
        }
        
        if (std::abs(psi_R - psi) < tolerance * psi) {
            R0 = R;
            break;
        }
    }
    
    if (Model::verbose && i >= limit) std::cout << "Inversion of the potential did not converge! Psi = " << psi << ", z = " << z << std::endl;
    return R0;
}

double Model::z_psi(double psi, double R, double tolerance, int limit) {
    if (psi == 0) return 1e306;
    
    double R2 = std::pow(R, 2);
    double z0 = 0, dz = 1.;
    int i = 0;
    while(i < limit) {
        i++;
        double z = z0 + dz;
        double z2 = std::pow(z, 2);
        double r = std::sqrt(R2 + z2);
        double psi_z = std::real(Model::psi(R2, z2, r));
        
        if (psi > psi_z) dz /= 2.;
        else {
            z0 = z;
            dz *= 2.;
        }
        
        if (std::abs(psi_z - psi) < tolerance * psi) {
            z0 = z;
            break;
        }
    }
    
    if (Model::verbose && i >= limit) std::cout << "Inversion of the potential did not converge! Psi = " << psi << ", R = " << R << std::endl;
    return z0;
}

