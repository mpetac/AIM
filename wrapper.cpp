#include <string>
#include <complex>
#include <iostream>
#include <fstream>
#include <time.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "src/Observables.hpp"
#include "src/Inversion.hpp"
#include "src/Model.hpp"
#include "src/Structs.hpp"
#include "src/baryons/Baryons.hpp"
#include "src/baryons/Baryons_H_2MN.hpp"
#include "src/halos/Halo.hpp"
#include "src/halos/Halo_NFW.hpp"
#include "src/halos/Halo_BUR.hpp" 
#include "src/halos/Halo_gNFW.hpp"
#include "src/halos/Halo_sABC.hpp"

//using namespace std;
namespace py = pybind11;

class PSDF {
    private:
        int verbose = 0;
        Halo_NFW *halo = NULL;
        Baryons *baryons = NULL;
        Model *model = NULL;
        Inversion *inversion = NULL;
        Observables *observables = NULL;
    public:
        PSDF(int v) {
            PSDF::verbose = v;
        }
        
        void setHalo(std::string type, py::array_t<double> p) {
            py::buffer_info params_buffer = p.request();
            double *params = (double *) params_buffer.ptr;
            
            if (type.compare("NFW") == 0) {
                halo_2p halo_params = {params[0], params[1]};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                Halo_NFW h(halo_params, halo_rot_params);
                PSDF::halo = &h;
            }
            /*
            else if (type.compare("BUR") == 0) {
                halo_2p halo_params = {params[0], params[1]};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                Halo_BUR h(halo_params, halo_rot_params);
                PSDF::halo = &h;
            }
            else if (type.compare("gNFW") == 0) {
                halo_6p halo_params = {params[0], params[1], 1., 3., params[4], 1.};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                Halo_gNFW h(halo_params, halo_rot_params);
                PSDF::halo = &h;
            }
            else if (type.compare("sABC") == 0) {
                halo_6p halo_params = {params[0], params[1], params[2], params[3], params[4], params[5]};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                Halo_sABC h(halo_params, halo_rot_params);
                PSDF::halo = &h;
            }
            */
            if (PSDF::verbose && PSDF::halo != NULL) std::cout << type << " halo successfully initialized!" << std::endl;
        }
        
        void setBaryons(std::string type, py::array_t<double> p) {
            py::buffer_info params_buffer = p.request();
            double *params = (double *) params_buffer.ptr;
            
            if (type.compare("H_2MN") == 0) {
                disk_3p disk1 = {params[0], params[1], params[2]};
                disk_3p disk2 = {params[3], params[4], params[5]};
                bulge_2p bulge = {params[6], params[7]};
                Baryons_H_2MN b(disk1, disk2, bulge);
                PSDF::baryons = &b;
            }
            if (PSDF::verbose && PSDF::baryons != NULL) std::cout << type << " baryons successfully initialized!" << std::endl;
        }
        
        void compute(int N, int M) {
            if (PSDF::halo == NULL) std::cout << "Halo not initialized!" << std::endl;
            else if (PSDF::baryons == NULL) std::cout << "Baryons not initialized!" << std::endl;
            else {
                std::cout << PSDF::halo << std::endl;
                std::cout << PSDF::halo->psi(0, 0, 0) << std::endl;
                
                /*
                Halo *h = (Halo *) halo;
                Baryons *b = (Baryons *) baryons;
                Model m(h, b);
                
                //std::cout << m.psi(0, 0, 0) << std::endl;
                
                
                Inversion inv(&m, N, M);
                Observables obs(&m, &inv);
                
                PSDF::model = &m;
                PSDF::inversion = &inv;
                PSDF::observables = &obs;
                */
                if (PSDF::verbose) std::cout << "PSDF computed!" << std::endl;
            }
        }
        
        double *rho(int N, py::array_t<double> R, py::array_t<double> z) {
            py::buffer_info Rpts_buffer = R.request();
            py::buffer_info zpts_buffer = z.request();
            double *Rpts = (double *) Rpts_buffer.ptr;
            double *zpts = (double *) zpts_buffer.ptr;
            double *result = new double[N];
            
            PSDF::observables->rho(N, Rpts, zpts, result);
            return result;
        }
        
        double v_mom(int mom, double R, double z) {
            return PSDF::observables->v_mom(mom, R, z);
        }
        
        double *pv_mag(int N, double R, double z) {
            double *result = new double[N];
            PSDF::observables->pv_mag(N, R, z, result);
            return result;
        }
        
        double *pv_merid(int N, double R, double z) {
            double *result = new double[N];
            PSDF::observables->pv_merid(N, R, z, result);
            return result;
        }
        
        double *pv_azim(int N, double R, double z) {
            double *result = new double[N];
            PSDF::observables->pv_azim(N, R, z, result);
            return result;
        }
        
        double *pv_rad(int N, double R, double z) {
            double *result = new double[N];
            PSDF::observables->pv_rad(N, R, z, result);
            return result;
        }
        
};


PYBIND11_MODULE(AIM, m) {
    py::class_<PSDF>(m, "PSDF")
        .def(py::init<bool>())
        .def("setHalo", &PSDF::setHalo)
        .def("setBaryons", &PSDF::setBaryons)
        .def("compute", &PSDF::compute)
        .def("rho", &PSDF::rho)
        .def("v_mom", &PSDF::v_mom)
        .def("pv_mag", &PSDF::pv_mag)
        .def("pv_merid", &PSDF::pv_merid)
        .def("pv_azim", &PSDF::pv_azim)
        .def("pv_rad", &PSDF::pv_rad);
}
