#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
namespace py = pybind11; 
using py::arg;

typedef py::array_t<double> dvec;

extern "C"{

void utm_geo_(double *rlon4,double *rlat4,double *rx4,double *ry4,
                const int *UTM_PROJECTION_ZONE,const int *iway);
}


auto utm_projection(const dvec &x,const dvec &y,int UTM_ZONE,const std::string &type = "GEO2UTM")
{
    assert(type == "GEO2UTM" || type == "UTM2GEO");
    int n = x.size();
    double xin[n],yin[n];
    memcpy(xin,x.data(),sizeof(double)*n);
    memcpy(yin,y.data(),sizeof(double)*n);
    dvec xout(n),yout(n);
    int iway;
    if(type == "UTM2GEO") {
        iway = 1;
        for(int i=0;i<n;i++){
            utm_geo_(xout.mutable_data()+i,yout.mutable_data()+i,
                        xin+i,yin+i,
                        &UTM_ZONE,&iway);
        }
        
    }
    else{
        iway = 0;
        for(int i=0;i<n;i++){
            utm_geo_(xin+i,yin+i,
                       xout.mutable_data()+i,yout.mutable_data()+i,
                        &UTM_ZONE,&iway);
        }
    }

    auto tp = std::make_tuple(xout,yout);

    return tp;
}

PYBIND11_MODULE(libutm,m){
    m.doc() = "UTM Projection Wrapper for Python";

    m.def("utm_projection",&utm_projection,arg("xin"),arg("yin"),
        arg("UTM_ZONE"),arg("projection_type") = "GEO2UTM",
        "UTM Projection Wrapper for Python");
}