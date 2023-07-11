#pragma once

class complex{
private:
    double x;
    double y;
public:
    complex(double x=0, double y=0);

    //Getters
    double Re() const{return x;}
    double Im() const{return y;}

    //Operations
    complex& operator+=(complex const&);
    complex& operator-=(complex const&);
    complex& operator*=(complex const&);
    complex operator*(complex const&) const;
    complex operator+(complex const&) const;
    complex operator-() const;
    complex operator-(complex const&) const;

    //Methodes
    complex conjugate() const;
    complex operator~() const;
};