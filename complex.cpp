#include "complex.h"

complex::complex(double x, double y) : x(x), y(y){}

bool complex::operator==(complex const& other) const{
    return x==other.x and y==other.y;
}

bool complex::operator!=(complex const& other) const{
    return not(*this==other);
}

complex& complex::operator*=(complex const& other){
    x = x*other.x - y*other.y;
    y = x*other.y + y*other.x;
    return *this;
}

complex& complex::operator+=(complex const& other){
    x+=other.x;
    y+=other.y;
    return *this;
}

complex& complex::operator-=(complex const& other){
    y-=other.y;
    x-=other.x;
}

complex complex::operator*(complex const& other) const{
    complex result(*this);
    return result*=other;
}

complex complex::operator+(complex const& other) const{
    complex result(*this);
    return result+=other;
}

complex complex::operator-() const{
    complex result;
    return result-=*this;
}

complex complex::operator-(complex const& other) const{
    complex result(*this);
    return result-=other;
}

complex complex::conjugate() const{
    return complex(x,-y);
}

complex complex::operator~() const{
    return this->conjugate();
}