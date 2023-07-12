#include "overloads.h"

std::string myTo_string(const std::complex<double>& c) {
    std::stringstream oss;

    double realPart(c.real());
    double imagPart = c.imag();

    if (realPart!=0 or imagPart==0) oss << c.real();

    if (imagPart!=0){
        if (imagPart >= 0) {
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
    return std::complex<double>(c.real(),-c.imag());
}