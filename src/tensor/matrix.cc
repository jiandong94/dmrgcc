#include "tensor/matrix.h"

// constructor
template<class T>
Matrix<T>::Matrix()
{
    row_ = 0;
    column_ = 0;
    total_element_num_ = 0;

    matrix_element_ = nullptr;
}

// constructor
//
template<class T>
Matrix<T>::Matrix(int row, int column)
{
    row_ = row;
    column_ = column;
    total_element_num_ = row*column;

    if(total_element_num_ > 0)
    {
        matrix_element_ = new T[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 0.0;
    }
    else matrix_element_ = nullptr;

}

// copy
template<class T>
Matrix<T>::Matrix(Matrix* tmp_matrix)
{
    row_ = tmp_matrix->row_;
    column_ = tmp_matrix->column_;
    total_element_num_ = tmp_matrix->total_element_num_;

    if(total_element_num_ > 0)
    {
        matrix_element_ = new T[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = tmp_matrix->matrix_element_[i];
    }
    else matrix_element_ = nullptr;

}


// destructor
//
template<class T>
Matrix<T>::~Matrix()
{
    delete[] matrix_element_;
}

// get row number
//
template<class T>
int Matrix<T>::get_row()
{
    return row_;
}

// get column number
//
template<class T>
int Matrix<T>::get_column()
{
    return column_;
}


// get total element number
//
template<class T>
int Matrix<T>::get_total_element_num()
{
    return total_element_num_;
}

// set matrix element
//
template<class T>
void Matrix<T>::set_matrix_element(int row, int column, T element)
{
    if(row >= row_ || column >= column_)
    {
        cout << "Matrix indices are out of range in set_matrix_element" << endl;
        exit(-1);
    }
    matrix_element_[row*column_+column] = element;
    
}

// get a matrix element
//
template<class T>
T Matrix<T>::get_matrix_element(int row, int column)
{
    if(row >= row_ || column >= column_)
    {
        cout << "Matrix indices are out of range in get_matrix_element" << endl;
        exit(-1);
    }

    return matrix_element_[row*column_+column];
}

// get matrix elements
//
template<class T>
T* Matrix<T>::get_matrix_element()
{
    return matrix_element_;
}

template<class T>
void Matrix<T>::PrintMatrix()
{
    if(total_element_num_!=0)
    {
        cout << "------------------------------"<< endl;
        cout << "Matrix Size: (" << row_ << "," << column_ << ")" << endl;
        for(int i=0;i<total_element_num_;++i)
        {
            if((i+1)%column_==0) cout << matrix_element_[i] << endl;
            else cout << matrix_element_[i] << ", ";
        }
        cout << endl;
    }
    else
    {
        cout << "The matrix_element_ is nullptr!" << endl;
    }
}

template<class T>
void Matrix<T>::WriteMatrix(const char* matrix_name)
{
    ofstream matrix_file;

    matrix_file.open(matrix_name, ios::binary|ios::out);

    WriteMatrix(matrix_file);

    matrix_file.close();
}

template<class T>
void Matrix<T>::WriteMatrix(ofstream &matrix_file)
{
    matrix_file.write((char*) &row_, sizeof(int));
    matrix_file.write((char*) &column_, sizeof(int));
    matrix_file.write((char*) &total_element_num_, sizeof(int));

    if(matrix_element_ != nullptr)
    {
        for(int i=0;i<total_element_num_;++i) 
            matrix_file.write((char*) &matrix_element_[i], sizeof(T));
    }
}

template<class T>
void Matrix<T>::ReadMatrix(const char* matrix_name)
{
    ifstream matrix_file;

    matrix_file.open(matrix_name, ios::binary|ios::in);

    ReadMatrix(matrix_file);

    matrix_file.close();
}

template<class T>
void Matrix<T>::ReadMatrix(ifstream &matrix_file)
{
    matrix_file.read((char*) &row_, sizeof(int));
    matrix_file.read((char*) &column_, sizeof(int));
    matrix_file.read((char*) &total_element_num_, sizeof(int));

    delete[] matrix_element_;
    if(total_element_num_ > 0)
    {
        matrix_element_ = new T[total_element_num_];
        for(int i=0;i<total_element_num_;++i)
            matrix_file.read((char*) &matrix_element_[i], sizeof(T));
    }
}


// set row and column to zero
//
template<class T>
void Matrix<T>::ResetMatrix()
{
    row_ = 0;
    column_ = 0;
    total_element_num_ = 0;

    matrix_element_ = nullptr;
}

// set matrix elements to zero
//
template<class T>
void Matrix<T>::ClearMatrix()
{
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 0.0;
}


// set matrix elements random
//
template<class T>
void Matrix<T>::RandomMatrix()
{
    for(int j=0;j<column_;++j) for(int i=0;i<row_;++i)
    {
        if(std::is_same<T, double>::value)
            set_matrix_element(i, j, 0.01*(rand()%100));
        else if(std::is_same<T, Complex>::value)
            set_matrix_element(i, j, Complex(0.01*(rand()%100), 0.01*(rand()%100)));
    }
    //for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 0.01*(rand()%100);

}

template<class T>
double Matrix<T>::SumSquareMatrix()
{
    double result = 0.0;
    for(int i=0;i<total_element_num_;++i)
    {
        result += SquareNorm(matrix_element_[i]);
    }
    return result;
}


template<class T>
void Matrix<T>::AddToMatrixElement(int row, int column, T element)
{
    matrix_element_[row*column_+column] += element;
}

template<class T>
void Matrix<T>::AddToMatrix(Matrix* tmp_matrix)
{
    if(tmp_matrix == nullptr || tmp_matrix->total_element_num_ == 0) 
    {
        cout << "tmp_matrix is nullptr in AddToMatrix" << endl;
        return;
    }
    if(total_element_num_ == 0)
    {
        row_ = tmp_matrix->row_;
        column_ = tmp_matrix->column_;
        total_element_num_ = tmp_matrix->total_element_num_;
        
        delete[] matrix_element_;
        matrix_element_ = new T[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = tmp_matrix->matrix_element_[i];
    }
    else
    {
        if(row_ != tmp_matrix->row_ || column_ != tmp_matrix->column_)
        {
            cout << "Matrix indeces do not match in AddToMatrix" << endl;
            exit(-1);
        }
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] += tmp_matrix->matrix_element_[i];

    }
}


template<class T>
void Matrix<T>::AddToMatrix(T factor, Matrix* tmp_matrix)
{
    if(tmp_matrix == nullptr || tmp_matrix->total_element_num_ == 0) 
    {
        cout << "tmp_matrix is nullptr in AddToMatrix" << endl;
        return;
    }
    if(total_element_num_ == 0)
    {
        row_ = tmp_matrix->row_;
        column_ = tmp_matrix->column_;
        total_element_num_ = tmp_matrix->total_element_num_;
        
        delete[] matrix_element_;
        matrix_element_ = new double[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 
                                              factor*tmp_matrix->matrix_element_[i];
    }
    else
    {
        if(row_ != tmp_matrix->row_ || column_ != tmp_matrix->column_)
        {
            cout << "Matrix indeces do not match in AddToMatrix" << endl;
            exit(-1);
        }
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] += 
                                              factor*tmp_matrix->matrix_element_[i];

    }
}

template<class T>
void Matrix<T>::MultiplyToScalar(T scalar)
{
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] *= scalar;
}

template<class T>
Matrix* Matrix<T>::MultiplyToMatrix(Matrix* tmp_matrix)
{
    if(column_ != tmp_matrix->row_)
    {
        cout << "Matrix indeces do not match in MultiplyToMatrix" << endl;
        exit(-1);
    }
    CBLAS_LAYOUT layout = CblasRowMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_TRANSPOSE transb = CblasNoTrans;
    MKL_INT m,n,k,lda,ldb,ldc;
    double alpha=1.0,beta=0.0;
    Matrix* result_matrix = new Matrix<T>(row_, tmp_matrix->column_);
    m = row_;
    k = column_;
    n = tmp_matrix->column_;
    lda = k;
    ldb = n;
    ldc = n;
    if(std::is_same<T, double>::value)
        cblas_dgemm(layout, transa, transb, m, n, k, alpha, matrix_element_, lda, 
                tmp_matrix->matrix_element_, ldb, beta, result_matrix->matrix_element_, ldc);
    else if(std::is_same<T, Complex>::value)
        cblas_zgemm(layout, transa, transb, m, n, k, alpha, matrix_element_, lda, 
                tmp_matrix->matrix_element_, ldb, beta, result_matrix->matrix_element_, ldc);
    
    return result_matrix;
}

template<class T>
Matrix* Matrix<T>::MatrixKronProduct(Matrix* tmp_matrix)
{
    Matrix* result_matrix;
    int row, column, position[2];
    double element[2];
    
    row = row_*tmp_matrix->get_row();
    column = column_*tmp_matrix->get_column();
    if(row==0 || column==0) return nullptr;
    result_matrix = new Matrix<T>(row, column);
    
    for(int i1=0;i1<row_;++i1) for(int j1=0;j1<column_;++j1)
    {
        element[0] = get_matrix_element(i1,j1);
        for(int i2=0;i2<tmp_matrix->get_row();++i2) for(int j2=0;j2<tmp_matrix->get_column();++j2)
        {
            position[0] = i2+i1*tmp_matrix->get_row();
            position[1] = j2+j1*tmp_matrix->get_column();
            element[1] = tmp_matrix->get_matrix_element(i2,j2);
            result_matrix->set_matrix_element(position[0], position[1], element[0]*element[1]);
        }

    }
    return result_matrix;
}

template<class T>
void Matrix<T>::MatrixElementProduct(Matrix* tmp_matrix)
{
    if(row_ != tmp_matrix->row_ || column_ != tmp_matrix->column_)
    {
        cout << "Matrix indeces do not match in MatrixElementProduct" << endl;
        exit(-1);
    }
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] *= tmp_matrix->matrix_element_[i];
}

template<class T>
void Matrix<T>::ChangeMatrix(int leigh, int truncate_dim)
{
    int tmp_total_element_num;
    T *tmp_matrix_element;
    if(leigh==0 && row_!=truncate_dim)
    {
        tmp_total_element_num = column_*truncate_dim;
        tmp_matrix_element = new T[tmp_total_element_num];
        for(int i=0;i<tmp_total_element_num;++i) 
            tmp_matrix_element[i] = 0.0;
        if(total_element_num_ <= tmp_total_element_num)
        {
            for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
                tmp_matrix_element[i*column_+j] = matrix_element_[i*column_+j];
        }
        row_ = truncate_dim;
        delete[] matrix_element_;
        total_element_num_ = tmp_total_element_num;
        matrix_element_ = tmp_matrix_element;

    }
    else if(leigh==1 && column_!=truncate_dim)
    {
        tmp_total_element_num = row_*truncate_dim;
        tmp_matrix_element = new T[tmp_total_element_num];
        for(int i=0;i<tmp_total_element_num;++i) 
            tmp_matrix_element[i] = 0.0;
        if(total_element_num_ <= tmp_total_element_num)
        {
            for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
                tmp_matrix_element[i*truncate_dim+j] = matrix_element_[i*column_+j];
        }
        column_ = truncate_dim;
        delete[] matrix_element_;
        total_element_num_ = tmp_total_element_num;
        matrix_element_ = tmp_matrix_element;
    }
    
}


template<class T>
Matrix* Matrix::TransposeMatrix()
{
    Matrix* tmp_matrix = new Matrix<T>(column_, row_);
    if(total_element_num_>0)
    {
        for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
            tmp_matrix->set_matrix_element(j, i, get_matrix_element(i, j));
    }
    return tmp_matrix;
}

template<class T>
Matrix* Matrix::HermitMatrix()
{
    Matrix* tmp_matrix = new Matrix<T>(column_, row_);
    if(total_element_num_>0)
    {
        for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
            tmp_matrix->set_matrix_element(j, i, Conj(get_matrix_element(i, j)));
    }
    return tmp_matrix;
}

template<class T>
Matrix* Matrix::ReshapeMatrix(int row, int column)
{
    RealMatrix* tmp_matrix;
    if(row*column != row_*column_)
    {
        cout << "Reshape dimension is not match in ReshapeMatrix" << endl;
        exit(-1);
    }
    if(total_element_num_>0)
    {
        tmp_matrix = new RealMatrix(row, column);
        for(int i=0;i<total_element_num_;++i)
            tmp_matrix->matrix_element_[i] = matrix_element_[i];
    }
    return tmp_matrix;
}

void RealMatrix::ReplaceMatrix(RealMatrix* tmp_matrix)
{
    delete[] matrix_element_;
    row_ = tmp_matrix->row_;
    column_ = tmp_matrix->column_;
    total_element_num_ = tmp_matrix->total_element_num_;
    
    if(total_element_num_ > 0)
    {
        matrix_element_ = new double[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = tmp_matrix->matrix_element_[i];
    }
    else matrix_element_ = nullptr;
}

void RealMatrix::ParallelMatrix(int leigh, int &num_unparallel, int* &position_unparallel, 
        RealMatrix* &transfer_tensor)
{
    bool zero;
    double element[2], prefactor;
    int info;
    
    if(leigh == 0)
    {
        position_unparallel = new int[column_];
        for(int i=0;i<column_;++i)
            position_unparallel[i] = -1;
        position_unparallel[0] = 0;

        transfer_tensor = new RealMatrix(column_, column_);
        transfer_tensor->set_matrix_element(0, 0, 1.0);
        num_unparallel = 1;

        for(int r=0;r<column_;++r)
        {
            for(int k=0;k<num_unparallel;++k)
            {
                info = num_unparallel;
                prefactor = 0.0;
                zero = true;
                
                for(int l=0;l<row_;++l)
                {
                    element[0] = get_matrix_element(l, r);
                    element[1] = get_matrix_element(l, position_unparallel[k]);

                    isParallelElement(element, position_unparallel[k], info, prefactor, zero);

                    if(info == -1) break;
                }

                if(info != -1)
                {
                    transfer_tensor->set_matrix_element(k, r, prefactor);
                    break;
                }
            }
            if(info == -1)
            {
                position_unparallel[num_unparallel] = r;
                transfer_tensor->set_matrix_element(num_unparallel, r, 1.0);
                num_unparallel++;
            }
        }
    }
    else if(leigh == 1)
    {
        position_unparallel = new int[row_];
        for(int i=0;i<row_;++i)
            position_unparallel[i] = -1;
        position_unparallel[0] = 0;

        transfer_tensor = new RealMatrix(row_, row_);
        transfer_tensor->set_matrix_element(0, 0, 1.0);
        num_unparallel = 1;

        for(int l=0;l<row_;++l)
        {
            for(int k=0;k<num_unparallel;++k)
            {
                info = num_unparallel;
                prefactor = 0.0;
                zero = true;
                
                for(int r=0;r<column_;++r)
                {
                    element[0] = get_matrix_element(l, r);
                    element[1] = get_matrix_element(position_unparallel[k], r);

                    isParallelElement(element, position_unparallel[k], info, prefactor, zero);

                    if(info == -1) break;
                }

                if(info != -1)
                {
                    transfer_tensor->set_matrix_element(l, k, prefactor);
                    break;
                }
            }
            if(info == -1)
            {
                position_unparallel[num_unparallel] = l;
                transfer_tensor->set_matrix_element(l, num_unparallel, 1.0);
                num_unparallel++;
            }
        }
    }
}

void RealMatrix::ExpanMatrix(int flag, RealMatrix* tmp_matrix)
{
    double* res_matrix_element;
    int res_row, res_column, res_total_element_num;
    if(flag==0 && column_==tmp_matrix->column_)
    {
        res_row = row_ + tmp_matrix->row_;
        res_column = column_;
        res_total_element_num = res_row*res_column;
        res_matrix_element = new double[res_total_element_num];
        for(int i=0;i<res_total_element_num;++i)
            res_matrix_element[i] = 0.0;
        for(int i=0;i<total_element_num_;++i)
            res_matrix_element[i] = matrix_element_[i];
        if(tmp_matrix->get_matrix_element() != nullptr)
        {
            for(int i=0;i<tmp_matrix->total_element_num_;++i)
                res_matrix_element[total_element_num_+i] = tmp_matrix->matrix_element_[i];
        }
    }
    else if(flag==1 && row_==tmp_matrix->row_)
    {
        res_row = row_;
        res_column = column_+tmp_matrix->column_;
        res_total_element_num = res_row*res_column;
        res_matrix_element = new double[res_total_element_num];
        for(int i=0;i<res_total_element_num;++i)
            res_matrix_element[i] = 0.0;
        for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
            res_matrix_element[i*res_column+j] = matrix_element_[i*column_+j];
        if(tmp_matrix->get_matrix_element() != nullptr)
        {
            for(int i=0;i<tmp_matrix->row_;++i) for(int j=0;j<tmp_matrix->column_;++j)
            {
                res_matrix_element[i*res_column+j+column_] = \
                tmp_matrix->matrix_element_[i*tmp_matrix->column_+j];
            }
        }
    }
    else
    {
        cout << "row or column dimension is not match in ExpanMatrix" << endl;
        exit(-1);
    }

    delete[] matrix_element_;
    row_ = res_row;
    column_ = res_column;
    total_element_num_ = res_total_element_num;
    matrix_element_ = res_matrix_element;
}


// svd 
// don't forget to delete left_matrix, right_matrix, sigular_value
void RealMatrix::SVDMatrix(RealMatrix* &left_matrix, RealMatrix* &right_matrix, 
                           double* &singular_value, int &singular_dim)
{
    int matrix_layout = LAPACK_ROW_MAJOR;
    int info;
    char jobu = 'S';
    char jobvt = 'S';
    lapack_int m, n, lda, ldu, ldvt;
    m = row_;
    n = column_;
    lda = n;
    ldu = min(m, n);
    ldvt = n;
    
    double superb[min(m,n)-1];

    singular_dim = min(row_, column_);

    left_matrix = new RealMatrix(m, ldu);
    right_matrix = new RealMatrix(ldvt, n);
    //left_matrix = new RealMatrix(m, m);
    //right_matrix = new RealMatrix(n, n);
    singular_value = new double[singular_dim];

    info = LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, matrix_element_, 
                          lda, singular_value, left_matrix->matrix_element_, 
                          ldu, right_matrix->matrix_element_, ldvt, superb);
    if(info>0)
    {
        cout << "The algorithm computing SVD failed to converge." << endl;
        exit(-1);
    }
}
*/
