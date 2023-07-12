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
    void print(std::ostream&) const;
    Matrix transpose() const;
    Matrix conjugate() const;
    Matrix operator~() const; //returns conjugate of transpose : M*
    bool is_hermitian() const;
};

template<typename T>
Matrix<T> operator*(T const&, Matrix<T> const&);

template<typename T>
std::ostream& operator<<(std::ostream&, Matrix<T> const&);

//============================================================

#pragma once

template<typename T>
Matrix<T>::Matrix(std::vector<std::vector<T>> matrix_){
    size_t maxe(0);
    for(auto& v:matrix_){
        maxe=std::max(maxe,v.size());
    }
    std::vector<std::vector<T>> result(matrix_.size(), std::vector<T> (maxe,0));
    for(size_t i(0); i<matrix_.size(); ++i){
        for(size_t j(0); j<matrix_[i].size(); ++j){
            result[i][j]=matrix_[i][j];
        }
    }
    matrix=result;
}

template<typename T>
Matrix<T>::Matrix(size_t m, size_t n, T value) : matrix(std::vector<std::vector<T>>(m, std::vector<T> (n,value))){}

template<typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T> const& other) const{
    if(matrix[0].size()!=other.matrix.size()) throw "Impossible operation";
    size_t m(matrix.size());
    size_t p(other.matrix.size());
    size_t n(other.matrix[0].size());
    std::vector<std::vector<T>> result(matrix.size(), std::vector<T> (other.matrix[0].size(),myZero()));
    for(size_t i(0); i<m; ++i){
        for(size_t j(0); j<n; ++j){
            for(size_t k(0); k<p; ++k){
                result[i][j] = result[i][j] + matrix[i][k]*other.matrix[k][j];
            }
        }
    }
    return Matrix<T>(result);
}

template<typename T>
bool Matrix<T>::operator==(Matrix<T> const& other) const{
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    if(other.matrix.size()!=m or other.matrix[0].size()!=n) return false;
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<m; ++j){
            if (not(matrix[i][j]==other.matrix[i][j])) return false;
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::operator!=(Matrix<T> const& other) const{
    return not(*this==other);
}

template<typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T> const& other) const{
    if(matrix.size()!=other.matrix.size() or matrix[0].size()!=other.matrix[0].size()) throw "Impossible operation";
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    std::vector<std::vector<T>> result(m, std::vector<T> (n,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<m; ++j){
            result[i][j]= matrix[i][j]+other.matrix[i][j];
        }
    }
    return Matrix<T>(result);
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const{
    std::vector<std::vector<T>> result(matrix.size(), std::vector<T> (matrix[0].size(),0));
    for(size_t i(0); i<matrix.size(); ++i){
        for(size_t j(0); j<matrix[0].size(); ++j){
            result[i][j]=-matrix[i][j];
        }
    }
    return Matrix<T>(result);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T> const& other) const{
    return *this + (-other);
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(Matrix<T> const& other){
    *this = *this + other;
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T> const& other){
    *this = (*this)*(other);
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(Matrix<T> const& other){
    *this = *this - other;
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::scale(T const& c){
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    std::vector<std::vector<T>> result(m, std::vector<T> (n,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<m; ++j){
            result[i][j]= c*matrix[i][j];
        }
    }
    return result;
}

template<typename T>
void Matrix<T>::print(std::ostream& out) const{
    int maxElementWidth = 0;

    // Find the maximum width of any element in the matrix
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            int elementWidth = myTo_string(element).length();
            if (elementWidth > maxElementWidth) {
                maxElementWidth = elementWidth;
            }
        }
    }

    // Print the matrix with centralized elements
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            out << std::setw(maxElementWidth) << std::left << myTo_string(element) << " ";
        }
        out << std::endl;
    }
}

template<typename T>
Matrix<T> Matrix<T>::transpose() const{
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    std::vector<std::vector<T>> result(n, std::vector<T> (m,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<m; ++j){
            result[i][j]= matrix[j][i];
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::conjugate() const{
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    std::vector<std::vector<T>> result(m, std::vector<T> (n,0));
    for(size_t i(0); i<n; ++i){
        for(size_t j(0); j<m; ++j){
            result[i][j]= ~matrix[i][j];
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator~() const{
    return transpose().conjugate();
} 

template<typename T>
bool Matrix<T>::is_hermitian() const{
    return *this==(~*this);
}

//External

template<typename T>
Matrix<T> operator*(T const& scalar, Matrix<T> const& m){
    return Matrix<T>(m).scale(scalar);
}

template<typename T>
std::ostream& operator<<(std::ostream& out, Matrix<T> const& M){
    M.print(out);
    return out;
}