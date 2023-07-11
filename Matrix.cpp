#include "Matrix.h"
using namespace std;

Matrix::Matrix(vector<vector<complex>> matrix_){
    size_t maxe(0);
    for(auto& v:matrix_){
        maxe=max(maxe,v.size());
    }
    vector<vector<complex>> result(matrix_.size(), vector<complex> (maxe,0));
    for(size_t i(0); i<matrix_.size(); ++i){
        for(size_t j(0); j<matrix_[i].size(); ++j){
            result[i][j]=matrix_[i][j];
        }
    }
    matrix=result;
}

Matrix::Matrix(size_t m, size_t n, complex value) : matrix(vector<vector<complex>>(m, vector<complex> (n,value))){}

Matrix Matrix::operator*(Matrix const& other) const{
    if(matrix[0].size()!=other.matrix.size()) throw "Impossible operation";
    size_t n(other.matrix.size());
    vector<vector<complex>> result(matrix.size(), vector<complex> (other.matrix[0].size(),0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<n; ++j){
            for(size_t k(0); k<n; ++k){
                result[i][j] += matrix[i][k]*matrix[k][j];
            }
        }
    }
    return Matrix(result);
}

bool Matrix::operator==(Matrix const& other){
    bool b(true);
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<n; ++j){
            if 
        }
    }

}

Matrix Matrix::operator+(Matrix const& other) const{
    if(matrix.size()!=other.matrix.size() or matrix[0].size()!=other.matrix[0].size()) throw "Impossible operation";
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    vector<vector<complex>> result(m, vector<complex> (n,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<n; ++j){
            result[i][j]= matrix[i][j]+other.matrix[i][j];
        }
    }
    return Matrix(result);
}

Matrix Matrix::operator-() const{
    vector<vector<complex>> result(matrix.size(), vector<complex> (matrix[0].size(),0));
    for(size_t i(0); i<matrix.size(); ++i){
        for(size_t j(0); j<matrix[0].size(); ++j){
            result[i][j]=-matrix[i][j];
        }
    }
    return Matrix(result);
}

Matrix Matrix::operator-(Matrix const& other) const{
    return *this + (-other);
}

Matrix& Matrix::operator+=(Matrix const& other){
    *this = *this + other;
    return *this;
}

Matrix& Matrix::operator*=(Matrix const& other){
    *this = (*this)*(other);
    return *this;
}

Matrix& Matrix::operator-=(Matrix const& other){
    *this = *this - other;
    return *this;
}

Matrix& Matrix::scale(complex const& c){
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    vector<vector<complex>> result(m, vector<complex> (n,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<n; ++j){
            result[i][j]= c*matrix[i][j];
        }
    }
}

Matrix Matrix::transpose() const{
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    vector<vector<complex>> result(n, vector<complex> (m,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<n; ++j){
            result[i][j]= matrix[j][i];
        }
    }
}

Matrix Matrix::conjugate() const{
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    vector<vector<complex>> result(n, vector<complex> (m,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<n; ++j){
            result[i][j]= ~matrix[j][i];
        }
    }
}

Matrix Matrix::operator~() const{
    return transpose().conjugate();
} 

bool Matrix::is_hermitian() const{
    return 
}

//External

Matrix operator*(complex const& scalar, Matrix const& m){
    return Matrix(m).scale(scalar);
}