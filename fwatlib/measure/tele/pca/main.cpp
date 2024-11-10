#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<iostream>

namespace py = pybind11;
using py::arg;

extern "C" {
void seis_pca(const float *stf_collect,int nrec,int npts,float *stf);
}

typedef py::array_t<float> farray;

farray 
PCA(farray &stf_collect) {
    auto buf = stf_collect.request();
    int nrec = buf.shape[0], npts = buf.shape[1];

    farray stf(npts);

    seis_pca(stf_collect.data(),nrec,npts,stf.mutable_data());

    return stf;
}

PYBIND11_MODULE(libpca,m){
    m.doc() = "PCA wrapper\n";
    m.def("PCA",&PCA,arg("stf_collect"),
          "PCA c++ wrapper");
}