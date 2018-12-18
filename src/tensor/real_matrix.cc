#include "tensor/real_matrix.h"

// constructor
RealMatrix::RealMatrix()
{
    row_ = 0;
    column_ = 0;
    total_element_num_ = 0;

    matrix_element_ = nullptr;
}

// constructor
//
RealMatrix::RealMatrix(int row, int column)
{
    row_ = row;
    column_ = column;
    total_element_num_ = row*column;

    if(total_element_num_ > 0)
    {
        matrix_element_ = new double[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 0.0;
    }
    else matrix_element_ = nullptr;

}

// copy
RealMatrix::RealMatrix(RealMatrix* tmp_matrix)
{
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


// destructor
//
RealMatrix::~RealMatrix()
{
    delete[] matrix_element_;
}

// get row number
//
int RealMatrix::get_row()
{
    return row_;
}

// get column number
//
int RealMatrix::get_column()
{
    return column_;
}

// get total element number
//
int RealMatrix::get_total_element_num()
{
    return total_element_num_;
}

// set matrix element
//
void RealMatrix::set_matrix_element(int row, int column, double element)
{
    if(row >= row_ || column >= column_)
    {
        cout << "Matrix indices are out of range!" << endl;
        exit(-1);
    }
    matrix_element_[row*column_+column] = element;
    
}

// get a matrix element
//
double RealMatrix::get_matrix_element(int row, int column)
{
    if(row >= row_ || column >= column_)
    {
        cout << "Matrix indices are out of range!" << endl;
        exit(-1);
    }

    return matrix_element_[row*column_+column];
}

// get matrix elements
//
double* RealMatrix::get_matrix_element()
{
    return matrix_element_;
}


void RealMatrix::PrintMatrix()
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
}

void RealMatrix::WriteMatrix(char* matrix_name)
{
    ofstream matrix_file;

    matrix_file.open(matrix_name, ios::binary|ios::out);

    WriteMatrix(matrix_file);

    matrix_file.close();
}

void RealMatrix::WriteMatrix(ofstream &matrix_file)
{
    matrix_file.write((char*) &row_, sizeof(int));
    matrix_file.write((char*) &column_, sizeof(int));
    matrix_file.write((char*) &total_element_num_, sizeof(int));

    if(matrix_element_ != nullptr)
    {
        for(int i=0;i<total_element_num_;++i) 
            matrix_file.write((char*) &matrix_element_[i], sizeof(double));
    }
}

void RealMatrix::ReadMatrix(char* matrix_name)
{
    ifstream matrix_file;

    matrix_file.open(matrix_name, ios::binary|ios::in);

    ReadMatrix(matrix_file);

    matrix_file.close();
}

void RealMatrix::ReadMatrix(ifstream &matrix_file)
{
    matrix_file.read((char*) &row_, sizeof(int));
    matrix_file.read((char*) &column_, sizeof(int));
    matrix_file.read((char*) &total_element_num_, sizeof(int));

    delete[] matrix_element_;
    if(total_element_num_ > 0)
    {
        matrix_element_ = new double[total_element_num_];
        for(int i=0;i<total_element_num_;++i)
            matrix_file.read((char*) &matrix_element_[i], sizeof(double));
    }
}


// set row and column to zero
//
void RealMatrix::ResetMatrix()
{
    row_ = 0;
    column_ = 0;
    total_element_num_ = 0;

    matrix_element_ = nullptr;
}

// set matrix elements to zero
//
void RealMatrix::ClearMatrix()
{
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 0.0;
}


// set matrix elements random
//
void RealMatrix::RandomMatrix()
{
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] = 0.01*(rand()%100);

}

//
//
void RealMatrix::AddToMatrixElement(int row, int column, double element)
{
    matrix_element_[row*column_+column] += element;
}

void RealMatrix::AddToMatrix(RealMatrix* tmp_matrix)
{
    if(tmp_matrix == nullptr || tmp_matrix->total_element_num_ == 0) return;
    if(total_element_num_ == 0)
    {
        row_ = tmp_matrix->row_;
        column_ = tmp_matrix->column_;
        total_element_num_ = tmp_matrix->total_element_num_;
        
        delete[] matrix_element_;
        matrix_element_ = new double[total_element_num_];
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] = tmp_matrix->matrix_element_[i];
    }
    else
    {
        if(row_ != tmp_matrix->row_ || column_ != tmp_matrix->column_)
        {
            cout << "Matrix indeces do not match!" << endl;
            exit(-1);
        }
        for(int i=0;i<total_element_num_;++i) matrix_element_[i] += tmp_matrix->matrix_element_[i];

    }
}

void RealMatrix::MultiplyToScalar(double scalar)
{
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] *= scalar;
}

RealMatrix* RealMatrix::MultiplyToMatrix(RealMatrix* tmp_matrix)
{
    if(column_ != tmp_matrix->row_)
    {
        cout << "Matrix indeces do not match in MultiplyToMatrix";
        exit(-1);
    }
    CBLAS_LAYOUT layout = CblasRowMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_TRANSPOSE transb = CblasNoTrans;
    MKL_INT m,n,k,lda,ldb,ldc;
    double alpha=1.0,beta=0.0;
    RealMatrix* result_matrix = new RealMatrix(row_, tmp_matrix->column_);
    m = row_;
    k = column_;
    n = tmp_matrix->column_;
    lda = k;
    ldb = n;
    ldc = n;

    cblas_dgemm(layout, transa, transb, m, n, k, alpha, matrix_element_, lda, 
                tmp_matrix->matrix_element_, ldb, beta, result_matrix->matrix_element_, ldc);
    
    return result_matrix;
}

void RealMatrix::MatrixElementProduct(RealMatrix* tmp_matrix)
{
    if(row_ != tmp_matrix->row_ || column_ != tmp_matrix->column_)
    {
        cout << "Matrix indeces do not match!" << endl;
        exit(-1);
    }
    for(int i=0;i<total_element_num_;++i) matrix_element_[i] *= tmp_matrix->matrix_element_[i];
}

RealMatrix* RealMatrix::TransposeMatrix()
{
    RealMatrix* tmp_matrix = new RealMatrix(column_, row_);
    if(total_element_num_>0)
    {
        for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
            tmp_matrix->set_matrix_element(j, i, get_matrix_element(i, j));
    }
    return tmp_matrix;
}

RealMatrix* RealMatrix::ReshapeMatrix(int row, int column)
{
    RealMatrix* tmp_matrix;
    if(row*column != row_*column_)
    {
        cout << "Reshape dimension is not match!" << endl;
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

void RealMatrix::CombineMatrix(int flag, RealMatrix* tmp_matrix)
{
    double* res_matrix_element;
    int res_row, res_column, res_total_element_num;
    if(flag==0 && column_==tmp_matrix->column_)
    {
        res_row = row_ + tmp_matrix->row_;
        res_column = column_;
        res_total_element_num = res_row*res_column;
        res_matrix_element = new double[res_total_element_num];
        for(int i=0;i<total_element_num_;++i)
            res_matrix_element[i] = matrix_element_[i];
        for(int i=0;i<tmp_matrix->total_element_num_;++i)
            res_matrix_element[total_element_num_+i] = tmp_matrix->matrix_element_[i];
    }
    else if(flag==1 && row_==tmp_matrix->row_)
    {
        res_row = row_;
        res_column = column_+tmp_matrix->column_;
        res_total_element_num = res_row*res_column;
        res_matrix_element = new double[res_total_element_num];
        for(int i=0;i<row_;++i) for(int j=0;j<column_;++j)
            res_matrix_element[i*res_column+j] = matrix_element_[i*column_+j];
        for(int i=0;i<tmp_matrix->row_;++i) for(int j=0;j<tmp_matrix->column_;++j)
        {
            res_matrix_element[i*res_column+j+column_] = \
            tmp_matrix->matrix_element_[i*tmp_matrix->column_+j];
        }
    }
    else
    {
        cout << "row or column dimension is not match" << endl;
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
                           double* &sigular_value, int &sigular_dim)
{
    int matrix_layout = LAPACK_ROW_MAJOR;
    int info;
    char jobu = 'A';
    char jobvt = 'A';
    lapack_int m, n, lda, ldu, ldvt;
    m = row_;
    n = column_;
    lda = n;
    ldu = m;
    ldvt = n;
    
    double superb[min(m,n)-1];

    sigular_dim = row_;
    if(row_ < column_) sigular_dim = column_;


    left_matrix = new RealMatrix(m, m);
    right_matrix = new RealMatrix(n, n);
    sigular_value = new double[sigular_dim];

    info = LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, matrix_element_, 
                          lda, sigular_value, left_matrix->matrix_element_, 
                          ldu, right_matrix->matrix_element_, ldvt, superb);
    if(info>0)
    {
        cout << "The algorithm computing SVD failed to converge." << endl;
        exit(-1);
    }
}
































