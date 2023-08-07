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

//Getters and Setters ===================================================

template<typename T>
Matrix<T> Matrix<T>::row(size_t i) const{
    if(i>=numRows()) throw std::invalid_argument("Row index out of range");
    return Matrix<T>({matrix[i]});
}

template<typename T>
Matrix<T> Matrix<T>::column(size_t i) const{
    if(i>=numCols()) throw std::invalid_argument("Column index out of range");
    Matrix<T> col({this->transpose().matrix[i]});
    return col.transpose();
}

template<typename T>
Matrix<T>& Matrix<T>::row(size_t index, Matrix const& M){
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
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::column(size_t index, Matrix const& N){
    std::vector<std::vector<T>> result;
    *this = this->transpose();
    Matrix<T> M(N.transpose());
    row(index, M);
    *this = this->transpose();
    return *this;
}

template<typename T>
T Matrix<T>::loc(size_t i, size_t j) const{
    return matrix[i][j];
}

template<typename T>
Matrix<T>& Matrix<T>::addRow(Matrix<T> const& M){
    size_t m(M.numRows());
    for(size_t i(0); i<m; ++i){
        matrix.push_back(M.matrix[i]);
    }
    *this = Matrix<T>(matrix); 
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::addColumn(Matrix<T> const& M){
    *this = this->transpose().addRow(M.transpose()).transpose();
    return *this;
}

//Operations ========================================================

template<typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T> const& other) const{
    if(matrix[0].size()!=other.matrix.size()) throw std::invalid_argument("Impossible operation");
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
    if(matrix.size()!=other.matrix.size() or matrix[0].size()!=other.matrix[0].size()) throw std::invalid_argument("Impossible operation");
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
            result[i][j]=myZero()-matrix[i][j];
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

//Methods =================================================================

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
Matrix<T> Matrix<T>::identity(size_t n){
    Matrix<T> result(n,n,myZero());
    for(size_t i(0); i<n; ++i){
        result.matrix[i][i]=myUnit();
    }
    return result;
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
Matrix<T> Matrix<T>::gauss_elimination(bool clear_null_lines, bool complete_gaussification, std::vector<size_t>* ref_points, Matrix<T>* Mirror) const{
    bool valid(Mirror!=nullptr); //Performs operations that were performed on M on Mirror matrix too, can yield inverse if M is full rank
    size_t m(numRows());
    size_t n(numCols());
    Matrix<T> result(*this);
    Matrix<T> temp(1,1,myZero());
    Matrix<T> tempMir(1,1,myZero());
    T scalar(myZero());
    size_t i(0);
    size_t start(0);
    std::vector<size_t> unused_rows(m);
    std::iota(unused_rows.begin(), unused_rows.end(), 0);
    for(size_t j(0); j<n; ++j ){
        i=m+1;
        for(size_t l(start); l<m; ++l){
            if(not(isZero(result.matrix[l][j]))){
                i=l;
                if(ref_points!=nullptr) ref_points->push_back(j);
                break;
            }
        }
        if(i<m){
            temp = result.row(i);
            if(valid) tempMir= Mirror->row(i);
            for(size_t k(i+1); k<m; ++k){
                temp = result.row(i);
                scalar = result.matrix[k][j]*myInverse(result.matrix[i][j]);
                result.row(k, result.row(k)-scalar*temp);
                if(valid) Mirror->row(k, Mirror->row(k)-scalar*tempMir);
            }
            if(complete_gaussification){
                for(size_t k(0); k<i; ++k){
                    scalar = result.matrix[k][j]*myInverse(result.matrix[i][j]);
                    result.row(k, result.row(k)-scalar*temp);
                    if(valid) Mirror->row(k, Mirror->row(k)-scalar*tempMir);
                }
                scalar = myInverse(result.matrix[i][j]);
                temp = scalar*temp;
                if(valid) tempMir = scalar*tempMir; 
                result.row(i,temp);
                if(valid) Mirror->row(i,tempMir);
            }
            result.row(i, result.row(start));
            result.row(start, temp);
            if(valid){
                Mirror->row(i, Mirror->row(start));
                Mirror->row(start, tempMir);
            }
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
Matrix<T> Matrix<T>::inverse() const{
    Matrix<T> Inverse(identity(numRows()));
    gauss_elimination(false,true,nullptr,&Inverse);
    return Inverse;
}

template<typename T>
Matrix<T> Matrix<T>::gram_schmidt(size_t start, size_t end, bool normalize) const{
    Matrix<T> result(*this);
    result=result.transpose();
    Matrix<T> temp1(1,1,myZero());
    Matrix<T> temp2(1,1,myZero());
    Matrix<T> temp(1,1,myZero());
    for(size_t i(start); i<end; ++i){
        temp = result.row(i);
        for(size_t j(0); j<i; ++j){
            temp1 = result.row(i)*(~(result.row(j)));
            temp2 = result.row(j)*(~(result.row(j)));
            temp = temp - (temp1.matrix[0][0]*myInverse(temp2.matrix[0][0]))*result.row(j);
        }
        result.row(i, temp);
        if(normalize){
            temp1 = result.row(i)*(~result.row(i));
            result.row(i, myInverse(mySqrt(temp1.matrix[0][0]))*result.row(i));
        }
    }
    return result.transpose();
}

template<typename T>
Matrix<T> Matrix<T>::gram_schmidt(bool normalize) const{
    return gram_schmidt(0, numCols(), normalize);
}

template<typename T>
Matrix<T> Matrix<T>::orthogonal_base(bool normalize) const{
    Matrix<T> result(this->transpose().gauss_elimination(true).transpose().gram_schmidt(true));
    return result;
}

template<typename T>
T Matrix<T>::determinant() const{
    if(numRows()!=numCols()) throw std::invalid_argument("Matrix must be square");
    size_t n(numRows());
    T det(myZero());
    Matrix<T> G(gauss_elimination());
    for(size_t i(0); i<n; ++i){
        if(i==0) det=G.matrix[i][i];
        else det = det*G.matrix[i][i];
    }
    return det;
}

template<typename T>
Matrix<T> Matrix<T>::kernel() const{
    std::vector<size_t> ref_points;
    Matrix<T> G(gauss_elimination(false, true, &ref_points));
    if(ref_points.size()==0) return identity(numCols());
    size_t n(numCols());
    size_t start;
    size_t end;
    Matrix<T> K(n,n-ref_points.size(),myZero());
    size_t index(0);
    for(size_t r(0); r<ref_points.size(); ++r){
        start = ref_points[r];
        if(r+1 == ref_points.size()) end=n;
        else end = ref_points[r+1];
        if(r==0){
            for(size_t i(0); i<ref_points[r]; ++i){
                K.column(index, G.column(i));
                K.matrix[i][index]= myUnit();
                index+=1;
            }
        }
        for(size_t i(start+1); i<end; ++i){
            K.column(index, -G.column(i));
            K.matrix[i][index]= myUnit();
            index+=1;
        }
    }
    return K;
}

template<typename T>
Matrix<T> Matrix<T>::image() const{
    return orthogonal_base(false);
}

template<typename T>
Matrix<T> Matrix<T>::complete_base(bool orthonormal) const{
    Matrix<T> K(transpose().kernel());
    Matrix<T> result(*this);
    result.addColumn(K);
    if(not(orthonormal)) return result;
    Matrix<T> O(result.orthogonal_base(true));
    for(size_t i(numRows()); i<result.numRows(); ++i){
        result.row(i, O.row(i));
    }
    return result;
}

template<typename T>
void Matrix<T>::QRD(Matrix<T>& Q, Matrix<T>& R) const{
    //Add something if M=0
    T norm(myZero());
    std::vector<std::vector<T>> Rmatrix;
    std::vector<T> coefs;
    Matrix<T> temp(1,1,myZero());
    Matrix<T> temp1(1,1,myZero());
    Matrix<T> temp2(1,1,myZero());
    T scalar(myZero());
    Matrix<T> tempCol(1,1,myZero());
    Matrix<T> nulvector(numRows(),1,myZero());
    bool firstime(true);
    size_t n(0);
    for(size_t i(0); i<numCols(); ++i){
        coefs={};
        tempCol = column(i);
        if(not firstime) n = Q.numCols();
        for(size_t j(0); j<n; ++j){
            temp1 = (column(i).transpose())*Q.column(j).conjugate();
            temp2 = (Q.column(j).transpose())*Q.column(j).conjugate();
            if(isZero(temp2.matrix[0][0])) scalar=myZero();
            else scalar = temp1.matrix[0][0]*myInverse(temp2.matrix[0][0]);
            tempCol = tempCol - scalar*Q.column(j);
            coefs.push_back(scalar);
        }
        temp = (~tempCol)*tempCol;
        norm=temp.matrix[0][0];
        if(not isZero(norm)){
            if(firstime){
                Q=tempCol;
                firstime=false;
            }
            else Q.addColumn(tempCol);
            coefs.push_back(myUnit());
            Rmatrix.push_back(coefs);
        }
        else{
            if(firstime){
                Q=nulvector;
                firstime=false;
            }
            Q.addColumn(nulvector);
            Rmatrix.push_back(coefs);
        }
    }
    Matrix<T> tempR(Rmatrix);
    R = tempR.transpose();
    //normalizing
    for(size_t k(0); k<Q.numCols(); ++k){
        temp = (~Q.column(k))*Q.column(k);
        norm = mySqrt(temp.matrix[0][0]);
        if(not(isZero(norm))){
            Q.column(k, myInverse(norm)*Q.column(k));
            R.row(k, norm*R.row(k));
        }
    }
    //squaring R
    while(R.numRows()<R.numCols()){
        R.addRow(nulvector.transpose());
    }
}

template<typename T>
Matrix<T> Matrix<T>::QR_algo() const{
    if(numCols()!=numRows()) throw std::invalid_argument("Can't find eigenvalues of non-square matrix");
    Matrix<T> Q(1,1,myZero());
    Matrix<T> R(1,1,myZero());
    Matrix<T> A(*this);
    QRD(Q,R);
    for(size_t i(0); i<50; ++i){
        A=R*Q;
        A.QRD(Q,R);
    }
    A=R*Q;
    return A;
}

template<typename T>
Matrix<T> Matrix<T>::diagonal_base(std::vector<T>* spectre) const{
    Matrix<T> D(QR_algo());
    Matrix<T> result(1,1,myZero());
    std::vector<std::vector<T>> map;
    bool dupe(false);
    T temp;
    for(size_t i(0); i<D.numCols(); ++i){
        dupe=false;
        temp = D.loc(i,i);
        for(size_t l(0); l<map.size(); ++l){
            if(isZero(map[l][0]-temp)){
                map[l].push_back(temp);
                dupe=true;
                break;
            }
        }
        if(not dupe){
            map.push_back({temp});
        }
    }
    size_t indexzero(0);
    for(size_t i(0); i<map.size(); ++i){
        if(not(isZero(map[i][0]))){
            D = *this - map[i][0]*identity(numCols()); //using D as temporary
            if(i==0){
                result = D.kernel();
            }
            else{
                result.addColumn(D.kernel());
            }
            if(spectre!=nullptr){
                for(size_t j(0); j<map[i].size();  ++j){
                    spectre->push_back(map[i][j]);
                }
            }
        }
        else indexzero=i;
    }
    //We add the eigen vectors associated to zero at the end.
    D = *this - map[indexzero][0]*identity(numCols());
    if(indexzero==0){
        result = D.kernel();
    }
    else{
        result.addColumn(D.kernel());
    }
    if(spectre!=nullptr){
        for(size_t j(0); j<map[indexzero].size();  ++j){
            spectre->push_back(map[indexzero][j]);
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::SVD(Matrix<T>& P, Matrix<T>& Q) const{
    std::vector<T> spectre;
    Matrix<T> D(numRows(), numCols(), myZero());
    Q = (~(*this)) * (*this);
    Q = Q.diagonal_base(&spectre).orthogonal_base(true);
    for(size_t i(0); i<spectre.size(); ++i){
        if(not isZero(spectre[i])){
            D.matrix[i][i]=mySqrt(spectre[i]);
            if(i==0){
                P = myInverse(mySqrt(spectre[i]))*((*this) * Q.column(i));
            }
            else{
                P.addColumn(myInverse(mySqrt(spectre[i]))*((*this) * Q.column(i)));
            }
        }
        else break;
    }
    P = P.complete_base(true);
    Q = ~Q;
    return D;
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