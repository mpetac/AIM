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
        bool verbose = 0;
        Halo *halo = NULL;
        Baryons *baryons = NULL;
        Model *model = NULL;
        Inversion *inversion = NULL;
        Observables *observables = NULL;
        
        py::array_t<double> py_array(int N, double *array) {
            auto result = py::array_t<double>(N);
            py::buffer_info result_buffer = result.request();
            double *result_ptr = (double *) result_buffer.ptr;
            for (int i = 0; i < N; i++) {
                result_ptr[i] = array[i];
            }
            return result;
        }
        
        double *c_array(py::array_t<double> array) {
            py::buffer_info array_buffer = array.request();
            double *result = (double *) array_buffer.ptr;
            return result;
        }
        
    public:
        PSDF(bool v) {
            PSDF::verbose = v;
        }
        
        void setHalo(std::string type, py::array_t<double> p) {
            double *params = PSDF::c_array(p);
            
            if (type.compare("NFW") == 0) {
                halo_2p halo_params = {params[0], params[1]};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                PSDF::halo = new Halo_NFW(halo_params, halo_rot_params);
            }
            else if (type.compare("BUR") == 0) {
                halo_2p halo_params = {params[0], params[1]};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                PSDF::halo = new Halo_BUR(halo_params, halo_rot_params);
            }
            else if (type.compare("gNFW") == 0) {
                halo_6p halo_params = {params[0], params[1], 1., 3., params[2], 1.};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                PSDF::halo = new Halo_gNFW(halo_params, halo_rot_params);
            }
            else if (type.compare("sABC") == 0) {
                halo_6p halo_params = {params[0], params[1], params[3], params[4], params[2], params[5]};
                halo_rot_2p halo_rot_params = {params[6], params[7]};
                PSDF::halo = new Halo_sABC(halo_params, halo_rot_params);
            }
            
            if (PSDF::verbose && PSDF::halo != NULL) std::cout << type << " halo successfully initialized!" << std::endl;
        }
        
        void setBaryons(std::string type, py::array_t<double> p) {
            double *params = PSDF::c_array(p);
            
            if (type.compare("H_2MN") == 0) {
                disk_3p disk1 = {params[0], params[1], params[2]};
                disk_3p disk2 = {params[3], params[4], params[5]};
                bulge_2p bulge = {params[6], params[7]};
                PSDF::baryons = new Baryons_H_2MN(disk1, disk2, bulge);
            }
            if (PSDF::verbose && PSDF::baryons != NULL) std::cout << type << " baryons successfully initialized!" << std::endl;
        }
        
        void compute(int N, int M) {
            if (PSDF::halo == NULL) std::cout << "Error! Halo not initialized!" << std::endl;
            else if (PSDF::baryons == NULL) std::cout << "Error! Baryons not initialized!" << std::endl;
            else {
                PSDF::model = new Model(PSDF::halo, PSDF::baryons);
                PSDF::inversion = new Inversion(PSDF::model, N, M, 1e-3, PSDF::verbose);
                PSDF::observables = new Observables(PSDF::model, PSDF::inversion, PSDF::verbose);
            }
        }
        
        py::array_t<double> rho(int N, py::array_t<double> R, py::array_t<double> z, double tolearance=1e-3) {
            py::buffer_info Rpts_buffer = R.request();
            py::buffer_info zpts_buffer = z.request();
            double *Rpts = (double *) Rpts_buffer.ptr;
            double *zpts = (double *) zpts_buffer.ptr;
            double *result = new double[N];

            PSDF::observables->rho(N, Rpts, zpts, result, tolearance);
            
            return PSDF::py_array(N, result);
        }
        
        py::array_t<double> rho_true(int N, py::array_t<double> R, py::array_t<double> z) {
            py::buffer_info Rpts_buffer = R.request();
            py::buffer_info zpts_buffer = z.request();
            double *Rpts = (double *) Rpts_buffer.ptr;
            double *zpts = (double *) zpts_buffer.ptr;
            double *result = new double[N];

            for (int i = 0; i < N; i++) {
                double R2 = std::pow(Rpts[i], 2), z2 = std::pow(zpts[i], 2);
                double r = std::sqrt(R2 + z2);
                result[i] = std::real(PSDF::model->rho(R2, z2, r));
            }
            
            return PSDF::py_array(N, result);
        }
        
        double v_mom(int mom, double R, double z, double tolearance=1e-3) {
            return PSDF::observables->v_mom(mom, R, z, tolearance);
        }
        
        py::array_t<double> pv_mag(int N, double R, double z, double tolearance=1e-3) {
            double *result = new double[2 * N];
            PSDF::observables->pv_mag(N, R, z, result, tolearance);
            return PSDF::py_array(2 * N, result);
        }
        
        py::array_t<double> pv_merid(int N, double R, double z, double tolearance=1e-3) {
            double *result = new double[2 * N];
            PSDF::observables->pv_merid(N, R, z, result, tolearance);
            return PSDF::py_array(2 * N, result);
        }
        
        py::array_t<double> pv_azim(int N, double R, double z, double tolearance=1e-3) {
            double *result = new double[2 * N];
            PSDF::observables->pv_azim(N, R, z, result, tolearance);
            return PSDF::py_array(2 * N, result);
        }
        
        py::array_t<double> pv_rad(int N, double R, double z, double tolearance=1e-3) {
            double *result = new double[2 * N];
            PSDF::observables->pv_rad(N, R, z, result, tolearance);
            return PSDF::py_array(2 * N, result);
        }
        
        py::array_t<double> pv_rel(int N, double R, double z, double tolearance=1e-1) {
            double *result = new double[2 * N];
            PSDF::observables->pv_rel(N, R, z, result, tolearance);
            return PSDF::py_array(2 * N, result);
        }
        
        py::array_t<double> occupation(int N_E, int N_Lz, py::array_t<double> Epts, py::array_t<double> Lzpts, double tolearance=1e-1) {
            int Npts = (N_E - 1) * (N_Lz - 1);
            py::buffer_info Epts_buffer = Epts.request();
            py::buffer_info Lzpts_buffer = Lzpts.request();
            double *E = (double *) Epts_buffer.ptr;
            double *Lz = (double *) Lzpts_buffer.ptr;
            double *result = new double[Npts];
            PSDF::observables->occupation(N_E, N_Lz, E, Lz, result, tolearance);
            return PSDF::py_array(Npts, result);
        }
        
        py::array_t<double> dd(int N, double t, double power, double R=8.122, double vR=11., double vPhi=-1., double vz=7., double vEarth=30., double vmax=1000., double tolearance=1e-3) {
            double *result = new double[N];
            PSDF::observables->dd(N, result, t, power, R, vR, vPhi, vz, vEarth, vmax, tolearance);
            return PSDF::py_array(N, result);
        }
        
};


PYBIND11_MODULE(AIM, m) {
    py::class_<PSDF>(m, "PSDF")
        .def(py::init<bool>())
        .def("setHalo", &PSDF::setHalo)
        .def("setBaryons", &PSDF::setBaryons)
        .def("compute", &PSDF::compute)
        .def("rho", &PSDF::rho)
        .def("rho_true", &PSDF::rho_true)
        .def("pv_mag", &PSDF::pv_mag)
        .def("pv_merid", &PSDF::pv_merid)
        .def("pv_azim", &PSDF::pv_azim)
        .def("pv_rad", &PSDF::pv_rad)
        .def("pv_rel", &PSDF::pv_rel)
        .def("v_mom", &PSDF::v_mom)
        .def("occupation", &PSDF::occupation)
        .def("dd", &PSDF::dd);
}
