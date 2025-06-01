#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <string>

namespace py = pybind11;
using py::arg;

const auto FCST = (py::array::c_style | py::array::forcecast) ;
typedef py::array_t<double,FCST> dvec;

#include <string>

extern "C" 
void measure_adj_fwat_(
    const double* obs_in,
    const double* syn_in,
    int npts,
    double t0,
    double dt,
    bool verbose,
    double tstart,
    double tend,
    double tt,
    double dtt,
    int nn,
    const char* chan_c,
    int imeas0_in,
    double TLONG_in,
    double TSHORT_in,
    bool RUN_BANDPASS_in,
    bool DISPLAY_DETAILS_in,
    bool OUTPUT_MEASUREMENT_FILES_in,
    bool COMPUTE_ADJOINT_SOURCE_in,
    double TSHIFT_MIN_in,
    double TSHIFT_MAX_in,
    double DLNA_MIN_in,
    double DLNA_MAX_in,
    double CC_MIN_in,
    int ERROR_TYPE_in,
    double DT_SIGMA_MIN_in,
    double DLNA_SIGMA_MIN_in,
    int ITAPER_in,
    double WTR_in,
    double NPI_in,
    double DT_FAC_in,
    double ERR_FAC_in,
    double DT_MAX_SCALE_in,
    double NCYCLE_IN_WINDOW_in,
    bool USE_PHYSICAL_DISPERSION_in,
    double* window_chi,
    double* tr_chi,
    double* am_chi,
    double* adj_out
);


auto measure(double t0,double dt, int npts,
            double tt,double dtt, int nn,
            double tstart,double tend, int imeas,
            double tlong,double tshort,bool verbose,
            const dvec &obs,const dvec &syn,
            bool RUN_BANDPASS,
            bool DISPLAY_DETAILS,
            bool OUTPUT_MEASUREMENT_FILES,
            bool COMPUTE_ADJOINT_SOURCE,
            double TSHIFT_MIN,
            double TSHIFT_MAX,
            double DLNA_MIN,
            double DLNA_MAX,
            double CC_MIN,
            int ERROR_TYPE,
            double DT_SIGMA_MIN,
            double DLNA_SIGMA_MIN,
            int ITAPER,
            double WTR,
            double NPI,
            double DT_FAC,
            double ERR_FAC,
            double DT_MAX_SCALE,
            double NCYCLE_IN_WINDOW,
            bool USE_PHYSICAL_DISPERSION)
{
    dvec window_chi,adj_src;
    adj_src.resize({nn}); window_chi.resize({20});
    double tr_chi,am_chi;
    char chan[4] = "BXZ";

    measure_adj_fwat(obs.data(),syn.data(),npts,t0,dt,verbose,tstart,tend,tt,dtt,nn,
                     chan,imeas,tlong,tshort,RUN_BANDPASS, DISPLAY_DETAILS, 
                     OUTPUT_MEASUREMENT_FILES,
                    COMPUTE_ADJOINT_SOURCE,
                    TSHIFT_MIN, TSHIFT_MAX, DLNA_MIN, DLNA_MAX, CC_MIN,
                    ERROR_TYPE, DT_SIGMA_MIN, DLNA_SIGMA_MIN,
                    ITAPER, WTR, NPI, DT_FAC,
                    ERR_FAC, DT_MAX_SCALE, NCYCLE_IN_WINDOW,
                    USE_PHYSICAL_DISPERSION,window_chi.mutable_data(),
                    &tr_chi,&am_chi,adj_src.mutable_data());

    return std::make_tuple(tr_chi,am_chi,window_chi,adj_src);
}

PYBIND11_MODULE(libmeas,m){
    m.doc() = "measure_adj python wrapper\n";
    m.def("measure",&measure,
            py::arg("t0"),
        py::arg("dt"),
        py::arg("npts"),
        py::arg("tt"),
        py::arg("dtt"),
        py::arg("nn"),
        py::arg("tstart"),
        py::arg("tend"),
        py::arg("imeas"),
        py::arg("tlong"),
        py::arg("tshort"),
        py::arg("verbose"),
        py::arg("obs"),
        py::arg("syn"),
        py::arg("RUN_BANDPASS"),
        py::arg("DISPLAY_DETAILS"),
        py::arg("OUTPUT_MEASUREMENT_FILES"),
        py::arg("COMPUTE_ADJOINT_SOURCE"),
        py::arg("TSHIFT_MIN"),
        py::arg("TSHIFT_MAX"),
        py::arg("DLNA_MIN"),
        py::arg("DLNA_MAX"),
        py::arg("CC_MIN"),
        py::arg("ERROR_TYPE"),
        py::arg("DT_SIGMA_MIN"),
        py::arg("DLNA_SIGMA_MIN"),
        py::arg("ITAPER"),
        py::arg("WTR"),
        py::arg("NPI"),
        py::arg("DT_FAC"),
        py::arg("ERR_FAC"),
        py::arg("DT_MAX_SCALE"),
        py::arg("NCYCLE_IN_WINDOW"),
        py::arg("USE_PHYSICAL_DISPERSION"),
        "Measure adjoint source based on input parameters.");

}