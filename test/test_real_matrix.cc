#include "tensor/real_matrix.h"

int main()
{
    cout << "=================================" << endl;
    cout << "         Test RealMatrix         " << endl;
    cout << "=================================" << endl;
    
    cout << "1. Constructor" << endl;
    // constructor
    RealMatrix* matrix = new RealMatrix();
    delete matrix;

    // constructor
    int m=2,n=2;
    matrix = new RealMatrix(m,m);
    for(int i=0;i<m;++i) for(int j=0;j<m;++j)
        matrix->set_matrix_element(i,j,i*m+j+1);
    // PrintMatrix
    matrix->PrintMatrix();


    
    cout << endl;
    cout << "2. Get" << endl;
    //get
    int row = matrix->get_row();
    int column = matrix->get_column();
    int total_element_num = matrix->get_total_element_num();
    int position_0_0 = matrix->get_matrix_element(0,0);
    double* matrix_element = matrix->get_matrix_element();
    matrix->PrintMatrix();
    cout << "row = " << row << endl;
    cout << "column = " << column << endl;
    cout << "total_element_num = " << total_element_num << endl;
    cout << "position_0_0 =  " << position_0_0 << endl;
    cout << "position_0_0 =  " << matrix_element[0] << endl;

    
    
    cout << endl;
    cout << "3. WriteMatrix ReadMatrix" << endl;
    // write read
    char const *char_matrix = "matrix.dat";
    matrix->PrintMatrix();
    cout << "Write matrix ..." << endl;
    matrix->WriteMatrix(char_matrix);
    cout << "Read matrix ..." << endl;
    RealMatrix* matrix_read = new RealMatrix();
    matrix_read->ReadMatrix(char_matrix);
    matrix_read->PrintMatrix();
    delete matrix_read;
    
    
    cout << endl;
    cout << "4. ResetMatrix ClearMatrix RandomMatrix SumSquareMatrix" << endl;

    cout << "Clear matrix:" << endl;
    matrix->ClearMatrix();
    matrix->PrintMatrix();

    cout << "Random matrix:" << endl;
    matrix->RandomMatrix();
    matrix->PrintMatrix();

    cout << "Sum of squared matrix elements" << endl;
    cout << "SumSquareMatrix: " << matrix->SumSquareMatrix() << endl;

    cout << "Reset matrix:" << endl;
    matrix->ResetMatrix();
    matrix->PrintMatrix();

    cout << endl;
    cout << "5. AddToMatrixElement AddToMatrix MultiplyToScalar " << endl;
    cout << "   MatrixElementProduct MultiplyToMatrix" << endl;
    // add element
    // add matrix
    delete matrix;
    matrix = new RealMatrix(m, n);
    for(int i=0;i<m;++i) for(int j=0;j<n;++j)
        matrix->set_matrix_element(i,j,i*m+j+1);
    matrix->PrintMatrix();

    cout << "AddToMatrixElement: matrix[0,0]+1" << endl;
    matrix->AddToMatrixElement(0,0,1);
    matrix->PrintMatrix();

    cout << "AddToMatrix: matrix+matrix" << endl;
    matrix->AddToMatrix(matrix);
    matrix->PrintMatrix();

    cout << "MultiplyToScalar: matrix*2" << endl;
    matrix->MultiplyToScalar(2);
    matrix->PrintMatrix();

    cout << "MatrixElementProduct: matrix.*matrix" << endl;
    matrix->MatrixElementProduct(matrix);
    matrix->PrintMatrix();

    cout << "MultiplyToMatrix: matrix*matrix" << endl;
    RealMatrix* matrix_multiply;
    matrix_multiply = new RealMatrix(matrix->MultiplyToMatrix(matrix));
    matrix_multiply->PrintMatrix();
    
    delete matrix_multiply;
    
    // test MKL 
    /*
    cout << endl;
    cout << "6. Test MKL dgemm speed" << endl;
    cout << "Test MKL speed: random_matrix[3000,3000]*random[3000,3000]" << endl;
    for(int i=0;i<10;++i)
    {
        RealMatrix* matrix_test = new RealMatrix(3000,3000);
        matrix_test->RandomMatrix();
        double start_time = GetWallTime();
        matrix_multiply = matrix_test->MultiplyToMatrix(matrix_test);
        double end_time = GetWallTime();
        cout << "[Iteration "<< i+1  << " Time :" << (end_time-start_time) << endl;
        
        delete matrix_test;
        delete matrix_multiply;
    }
    */
    


    cout << endl;
    cout << "7. ChangeMatrix TransposeMatrix ReshapeMatrix ExpanMatrix" << endl;
    matrix->PrintMatrix();
    
    cout << "Transpose matrix: " << endl;
    RealMatrix* matrix_transpose = matrix->TransposeMatrix();
    matrix_transpose->PrintMatrix();

    cout << "Reshape matrix: (2,2)->(1,4)" << endl;
    RealMatrix* matrix_reshape = matrix->ReshapeMatrix(1,4);
    matrix_reshape->PrintMatrix();

    cout << "Change matrix: (2,2)->(2,3)" << endl;
    matrix->ChangeMatrix(1,3);
    matrix->PrintMatrix();

    cout << "Expan matrix: (2,3),(2,3)->(2,6)" << endl;
    matrix->ExpanMatrix(1,matrix);
    matrix->PrintMatrix();

    delete matrix_transpose, matrix_reshape;


    cout << endl;
    cout << "8. SVDMatrix" << endl;
    delete matrix;
    matrix = new RealMatrix(2,3);
    for(int i=0;i<2;++i) for(int j=0;j<3;++j)
        matrix->set_matrix_element(i,j,i*3+j+1);
    matrix->PrintMatrix();
    RealMatrix *left_matrix, *right_matrix;
    double* sigular_value;
    int sigular_dim;
    matrix->SVDMatrix(left_matrix, right_matrix, sigular_value, sigular_dim);
    left_matrix->PrintMatrix();
    right_matrix->PrintMatrix();
    cout << "sigular dim: " << sigular_dim << endl;
    cout << "sigular value: ";
    for(int i=0;i<sigular_dim;++i) 
        cout << sigular_value[i] << " ";
    cout << endl;

    
    delete left_matrix, right_matrix;
    delete[] sigular_value;
    delete matrix;


    return 0;
}
