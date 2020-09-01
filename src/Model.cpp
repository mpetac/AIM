#include "Model.hpp"

std::complex<double> Model::psi(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi(R2, z2, r) + Model::baryons->psi(R2, z2, r);
}

std::complex<double> Model::psi_dR2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_dR2(R2, z2, r) + Model::baryons->psi_dR2(R2, z2, r);
}

std::complex<double> Model::psi_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_dz2(R2, z2, r) + Model::baryons->psi_dz2(R2, z2, r);
}

std::complex<double> Model::psi_d2R2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_d2R2z2(R2, z2, r) + Model::baryons->psi_d2R2z2(R2, z2, r);
}

std::complex<double> Model::psi_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->psi_d2z2(R2, z2, r) + Model::baryons->psi_d2z2(R2, z2, r);
}

std::complex<double> Model::rho(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->rho(R2, z2, r);
}

std::complex<double> Model::rho_dz2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->rho_dz2(R2, z2, r);
}

std::complex<double> Model::rho_d2z2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    return Model::halo->rho_d2z2(R2, z2, r);
}

std::complex<double> Model::rho_d2psi2(std::complex<double> R2, std::complex<double> z2, std::complex<double> r) {
    std::complex<double> dpsi_dz2 = Model::psi_dz2(R2, z2, r);
    return Model::psi(R2, z2, r) * std::pow(dpsi_dz2, -2) - Model::rho_dz2(R2, z2, r) * Model::psi_d2z2(R2, z2, r) * std::pow(dpsi_dz2, -3);
}

std::complex<double> Model::v_phi(std::complex<double> R2) {
    return Model::halo->v_phi(R2);
}



bool Model::is_rotating() {
    return Model::halo->is_rotating();
}

double psi_inverse_diff(const gsl_vector* v, void* p) {
    psi_inverse_params *params = (psi_inverse_params *) p;
    Model *model = (Model *) params->model;
    std::complex<double> z2(gsl_vector_get(v, 0), gsl_vector_get(v, 1));
    std::complex<double> r = std::sqrt(params->R2 + z2);
    std::complex<double> psi_val = model->psi(params->R2, z2, r);
    return std::abs((params->xi - psi_val) / psi_val);
}

std::complex<double> Model::psi_inverse(std::complex<double> xi, double E, double Lz, std::complex<double> z0, double tolerance, int limit) {
    int i = 0;
    double diff = 1;
    if (std::abs(z0) > 1e4) limit = 1000;
    
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
    
    return result;
}


double Model::Rcirc(double E, double tolerance, int limit) {
    int i = 0;
    double Rc = 0, dR = 1;
    while(i < limit) {
        i++;
        double R = Rc + dR;
        double R2 = std::pow(R, 2);
        double E_R = std::real(Model::psi(R2, 0, R) + R2 * Model::psi_dR2(R2, 0, R));
        
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
    
    if (Model::verbose && i >= limit) std::cout << "Rc inversion did not converge! E = " << E << std::endl;
    return Rc;
}


