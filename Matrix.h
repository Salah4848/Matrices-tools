#pragma once

#include <complex>
#include <vector>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <algorithm>

template<typename T>
class Matrix{ //We will call the matrix of the current instance M
private:
    std::vector<std::vector<T>> matrix;
public:
    Matrix(std::vector<std::vector<T>> matrix);
    Matrix(std::initializer_list<std::initializer_list<T>> matrix);
    Matrix(std::size_t m, std::size_t n, T value);

    //Getters and Setters
    size_t numRows() const{return matrix.size();}
    size_t numCols() const{return matrix[0].size();}
    Matrix row(size_t index) const;
    Matrix column(size_t index) const;
    Matrix& row(size_t index, Matrix const&); //Replaces the indexed row with input matrix (extends matrix with zeros if size differs)
    Matrix& column(size_t index, Matrix const&); //Replaces the indexed column with input matrix if possible (extends matrix with zeros if size differs)
    T loc(size_t i, size_t j) const;
    Matrix& addRow(Matrix const&);
    Matrix& addColumn(Matrix const&);

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
    void print(std::ostream&) const; //Works only if myTo_string(T) is defined
    static Matrix identity(size_t n); // returns identity matrix of size numRows()
    Matrix transpose() const;
    Matrix conjugate() const;
    Matrix operator~() const; //returns conjugate of transpose : M*
    bool is_hermitian() const;
    //All the following methods need type T to define a field
    Matrix gauss_elimination(bool clear_null_lines=false, bool complete_gaussification = false, std::vector<size_t>* ref_points = nullptr, Matrix* Mirror=nullptr) const; //Uses gauss elemination on M. need myInverse defined defined.
    Matrix inverse() const; //Returns inverse of M if possible using gauss_eliminiation
    //All the following methods need type T to define the complex or real fields
    Matrix gram_schmidt(size_t start, size_t end, bool normalize=true) const; //Directly performs gram-schmidt algorithm on the columns from start index to end index
    Matrix gram_schmidt(bool normalize=true) const; //Takes start as 0 and end as numCols.
    Matrix orthogonal_base(bool normalize=true) const; //Uses Gram-Schmidt to find orthonormal base for the columns of M, works only if gauss_elimination works, normalization only if MySqrt is defined
    T determinant() const; //Requires gauss elimination
    Matrix kernel() const; //Returns matrix with columns that form a base for the kernel of M
    Matrix image() const; //Returns matrix forming base of Im(M)
    Matrix complete_base(bool orthonormal=false) const; //Adds columns to M such that rank(M)=numcols and new vectors are orthogonal to Im(M) (if orthonormal=true)
    void QRD(Matrix& Q, Matrix& R) const; //Yields QR decomp of M, where Q is unitary and R upper triangular
    Matrix QR_algo() const; //Applies QR algorithm on M and returns final matrix
    Matrix diagonal_base(std::vector<T>* spectre=nullptr) const; //returns diagonal base for Im(M)
    Matrix SVD(Matrix& P, Matrix& Q) const; //returns the diagonal matrix D of singular values of M. Modifies P and Q such that : M=PDQ. Works only if orthonormal_base() works.

};

//External functions

template<typename T>
Matrix<T> operator*(T const&, Matrix<T> const&);

template<typename T>
std::ostream& operator<<(std::ostream&, Matrix<T> const&);

template<typename T>
Matrix<std::complex<T>> operator*(T const& scalar, Matrix<std::complex<T>> const& matrix);

#include "Matrix.tpp"