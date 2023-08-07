#include "overloads.h"
#include "Matrix.h"
#include <iostream>
using namespace std;

int main(){
    complex<double> c1(1,1);
    complex<double> c2(1,-1);
    Matrix<complex<double>> M({{c1,2,4},{2.*c1,4,8},{1,0,0}});
    Matrix<complex<double>> P(1,1,myZero());
    Matrix<complex<double>> Q(1,1,myZero());
    Matrix<complex<double>> D(M.SVD(P,Q));
    cout<<(~P)*M*(~Q)<<"\n";
    cout<<D<<"\n";
    cout<<P;
}