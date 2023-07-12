#include "overloads.h"
#include "Matrix.h"
#include <iostream>
using namespace std;

int main(){
    complex<double> c1(1,1);
    complex<double> c2(1,-1);
    Matrix<complex<double>> M({{c1,0,1},{0,c2,0},{0,0,1}});
    cout<<M;
    cout<<"\n";
    cout<<~M;
    cout<<"\n";
    cout<<(~M)*M;

    return 0;
}