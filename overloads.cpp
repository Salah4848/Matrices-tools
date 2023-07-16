#include "overloads.h"

std::string myTo_string(const std::complex<double>& c) {
    std::stringstream oss;

    double realPart(c.real());
    double imagPart = c.imag();

    if (realPart!=0 or imagPart==0) oss << c.real();

    if (imagPart!=0){
        if (imagPart >= 0 and realPart!=0) {
            oss << "+";
        }

        oss << imagPart << "i";
    }
    return oss.str();
}

std::complex<double> myZero(){
    return std::complex<double>(0,0);
}

std::complex<double> operator~(std::complex<double> const& c){
    return std::complex<double>(c.real(),-(c.imag()));
}

std::complex<double> myInverse(std::complex<double> const& c){
    if(isZero(c)) throw "Division by zero";
    std::complex<double> one(1,0);
    return one/c;
}

bool isZero(std::complex<double> const& c){
    if(c.real()*c.real()+c.imag()*c.imag()<1e-12) return true;
    return false;
}

std::complex<double> mySqrt(std::complex<double> const& c){
    return std::complex<double> (std::sqrt(std::abs(c)), 0);
}

std::complex<double> myUnit(){
    return std::complex<double>(1,0);
}
