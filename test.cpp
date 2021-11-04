//#include<xtensor/xtensor.hpp>
//#include<xtensor/xarray.hpp>
#include<xtensor.hpp>
#include<xtensor-blas/xblas.hpp>
#include<xtensor-blas/xlinalg.hpp>
#include<iostream>
#include<time.h>

int main(){
    unsigned int n = 2000;
    xt::xtensor<double,2> a({n,n}),b({n,n}),c({n,n});
    //a = xt::zeros<double>({n,n});
    clock_t tic,toc;

    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        a(i,j) = 1.0 * i;
        b(i,j) = 1.0 * j;
    }}

    tic = clock();
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        double tmp = 0.0;
        for(int k=0;k<n;k++) tmp += a(i,k) * b(k,j);
        c(i,j) = tmp;
    }}
    toc = clock();
    printf("%f\n",(double)((toc - tic) / CLOCKS_PER_SEC));

    tic = clock();
    c = xt::linalg::dot(a,b);
    toc = clock();
    printf("%f\n",(double)((toc - tic) / CLOCKS_PER_SEC));
}