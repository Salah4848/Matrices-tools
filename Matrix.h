#pragma once

#include <vector>
#include "complex.h"

class Matrix{
private:
    std::vector<std::vector<complex>> matrix;
public:
    Matrix(std::vector<std::vector<complex>> matrix);
    Matrix(std::size_t m, std::size_t n, complex value=0);

    //Getters
    size_t numRows() const{return matrix.size();}
    size_t numCols() const{return matrix[0].size();}

    //operations
    bool operator==(Matrix const&) const;
    bool operator!=(Matrix const&) const;
    Matrix operator*(Matrix const&) const;
    Matrix operator+(Matrix const&) const;
    Matrix operator-() const;
    Matrix operator-(Matrix const&) const;
    Matrix& operator+=(Matrix const&);
    Matrix& operator*=(Matrix const&);
    Matrix& operator-=(Matrix const&);
    Matrix& scale(complex const&);

    //Methods
    Matrix transpose() const;
    Matrix conjugate() const;
    Matrix operator~() const; //returns conjugate of transpose : M*
    bool is_hermitian() const;
};

Matrix operator*(complex const&, Matrix const&);