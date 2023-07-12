#pragma once

#include <string>
#include <sstream>
#include <complex>

std::string myTo_string(std::complex<double> const& c);
std::complex<double> myZero();
std::complex<double> operator~(std::complex<double> const&);