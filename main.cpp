#include "overloads.h"
#include "Matrix.h"
#include <iostream>
using namespace std;

int main(){
    complex<double> c1(1,1);
    complex<double> c2(1,-1);
    Matrix<complex<double>> M({{c1,2,1},{4,c2,7},{1,0,1}});
    cout<<M<<"\n";
    cout<<M.gauss_elimination();
    return 0;
}