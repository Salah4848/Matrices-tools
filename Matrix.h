#pragma once

#include <vector>
#include <iomanip>
#include <iostream>

template<typename T>
class Matrix{
private:
    std::vector<std::vector<T>> matrix;
public:
    Matrix(std::vector<std::vector<T>> matrix);
    Matrix(std::size_t m, std::size_t n, T value);

    //Getters
    size_t numRows() const{return matrix.size();}
    size_t numCols() const{return matrix[0].size();}
    Matrix row(size_t index) const;
    Matrix column(size_t index) const;

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
    Matrix scale(T const&);

    //Methods
    void print(std::ostream&) const; //Wroks only if myTo_string(T) is defined
    Matrix transpose() const;
    Matrix conjugate() const;
    Matrix operator~() const; //returns conjugate of transpose : M*
    bool is_hermitian() const;
    Matrix orthonormal_base() const; //Uses Gram-Schmidt to find orthonormal base for the columns of M, works only if mySqrt(T) is defined
    Matrix SVD(Matrix& P, Matrix& Q) const; //returns the diagonal matrix D of singular values of M. Modifies P and Q such that : M=PDQ. Works only if orthonormal_base() works.

};

//External functions

template<typename T>
Matrix<T> operator*(T const&, Matrix<T> const&);

template<typename T>
std::ostream& operator<<(std::ostream&, Matrix<T> const&);

#include "Matrix.tpp"