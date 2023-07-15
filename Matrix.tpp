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
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> matrix_){
    matrix.reserve(matrix_.size());
    for (const auto& row : matrix_) {
        matrix.emplace_back(row.begin(), row.end());
    }
    *this = Matrix<T> (matrix);
}

template<typename T>
Matrix<T>::Matrix(size_t m, size_t n, T value) : matrix(std::vector<std::vector<T>>(m, std::vector<T> (n,value))){}

template<typename T>
Matrix<T> Matrix<T>::row(size_t i) const{
    return Matrix<T>({matrix[i]});
}

template<typename T>
Matrix<T> Matrix<T>::column(size_t i) const{
    Matrix<T> col({this->transpose().matrix[i]});
    return col.transpose();
}

template<typename T>
void Matrix<T>::row(size_t index, Matrix const& M){
    std::vector<std::vector<T>> result;
    size_t n(numRows());
    for(size_t i(0); i<n; ++i){
        if(i==index){
            for(auto const& r : M.matrix){
                result.push_back(r);
            }
        }
        else result.push_back(matrix[i]);
        }
    Matrix<T> temp(result);
    *this = temp;
}

template<typename T>
void Matrix<T>::column(size_t index, Matrix const& N){
    std::vector<std::vector<T>> result;
    *this = this->transpose();
    Matrix<T> M(N.transpose());
    row(index, M);
    *this = this->transpose();
}

template<typename T>
T Matrix<T>::loc(size_t i, size_t j) const{
    return matrix[i][j];
}

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
    return result;
}

template<typename T>
bool Matrix<T>::operator==(Matrix<T> const& other) const{
    size_t m(matrix.size());
    size_t n(matrix[0].size());
    if(other.matrix.size()!=m or other.matrix[0].size()!=n) return false;
    for(size_t i(0); i<m; ++i){
        for(size_t j(0); j<n; ++j){
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
    std::vector<std::vector<T>> result(m, std::vector<T> (n,myZero()));
    for(size_t i(0); i<m; ++i){
        for(size_t j(0); j<n; ++j){
            result[i][j]= matrix[i][j]+other.matrix[i][j];
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const{
    std::vector<std::vector<T>> result(matrix.size(), std::vector<T> (matrix[0].size(),myZero()));
    for(size_t i(0); i<matrix.size(); ++i){
        for(size_t j(0); j<matrix[0].size(); ++j){
            result[i][j]=-matrix[i][j];
        }
    }
    return result;
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
    std::vector<std::vector<T>> result(m, std::vector<T> (n,myZero()));
    for(size_t i(0); i<m; ++i){
        for(size_t j(0); j<n; ++j){
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
    std::vector<std::vector<T>> result(n, std::vector<T> (m,myZero()));
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
    std::vector<std::vector<T>> result(m, std::vector<T> (n,myZero()));
    for(size_t i(0); i<m; ++i){
        for(size_t j(0); j<n; ++j){
            result[i][j]=~matrix[i][j];
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

template<typename T>
Matrix<T> Matrix<T>::gauss_elimination(bool clear_null_lines) const{
    size_t m(numRows());
    size_t n(numCols());
    Matrix<T> result(*this);
    Matrix<T> temp(1,1,myZero());
    size_t i(0);
    size_t start(0);
    std::vector<size_t> unused_rows(m);
    std::iota(unused_rows.begin(), unused_rows.end(), 0);
    for(size_t j(0); j<n; ++j ){
        i=m+1;
        for(size_t l(start); l<m; ++l){
            if(not(isZero(result.matrix[l][j]))){
                i=l;
                break;
            }
        }
        if(i<m){
            temp = result.row(i);
            for(size_t k(i+1); k<m; ++k){
                temp = result.row(i);
                result.row(k, result.row(k)-(result.matrix[k][j]*myInverse(result.matrix[i][j]))*temp);
            }
            result.row(i, result.row(start));
            result.row(start, temp);
            start++;
        }
    }
    if(clear_null_lines){
        std::vector<std::vector<T>> clean(result.matrix.begin(), result.matrix.begin() + start);
        return clean;
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::orthogonal_base(bool normalize) const{
    Matrix<T> temp1(1,1,myZero());
    Matrix<T> temp2(1,1,myZero());
    Matrix<T> result(this->transpose().gauss_elimination(true));
    for(size_t i(0); i<result.numRows(); ++i){
        for(size_t j(0); j<i; ++j){
            temp1 = result.row(i)*(~(result.row(j)));
            temp2 = result.row(j)*(~(result.row(j)));
            result.row(i, result.row(i)-(temp1.matrix[0][0]*myInverse(temp2.matrix[0][0]))*result.row(j));
        }
        if(normalize){
            temp1 = result.row(i)*(~result.row(i));
            result.row(i, myInverse(mySqrt(temp1.matrix[0][0]))*result.row(i));
        }
    }
    return result.transpose();
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

template<typename T>
Matrix<std::complex<T>> operator*(T const& scalar, Matrix<std::complex<T>> const& matrix) {
    std::complex<T> complexScalar(scalar, 0.0);
    return complexScalar * matrix;
}