#include "tensor/complex_matrix.h"

int main()
{
    cout << "=================================" << endl;
    cout << "         Test ComplexMatrix         " << endl;
    cout << "=================================" << endl;
    
    cout << "1. Constructor" << endl;
    // constructor
    ComplexMatrix* matrix = new ComplexMatrix();
    delete matrix;

    // constructor
    int m=2,n=2;
    matrix = new ComplexMatrix(m,m);
    for(int i=0;i<m;++i) for(int j=0;j<m;++j)
        matrix->set_matrix_element(i,j,i*m+j+1);
        //matrix->set_matrix_element(i,j,Complex(i*m+j+1,i*m+j+1));
    // PrintMatrix
    matrix->PrintMatrix();


    
    cout << endl;
    cout << "2. Get" << endl;
    //get
    int row = matrix->get_row();
    int column = matrix->get_column();
    int total_element_num = matrix->get_total_element_num();
    Complex position_0_0 = matrix->get_matrix_element(0,0);
    Complex* matrix_element = matrix->get_matrix_element();
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
    ComplexMatrix* matrix_read = new ComplexMatrix();
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
    matrix = new ComplexMatrix(m, n);
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
    ComplexMatrix* matrix_multiply;
    matrix_multiply = new ComplexMatrix(matrix->MultiplyToMatrix(matrix));
    matrix_multiply->PrintMatrix();
    
    delete matrix_multiply;
    
    // test MKL 
    /* 
    cout << endl;
    cout << "6. Test MKL dgemm speed" << endl;
    cout << "Test MKL speed: random_matrix[3000,3000]*random[3000,3000]" << endl;
    for(int i=0;i<10;++i)
    {
        ComplexMatrix* matrix_test = new ComplexMatrix(3000,3000);
        matrix_test->RandomMatrix();
        double start_time = GetWallTime();
        //matrix_multiply = matrix_test->MultiplyToMatrix(matrix_test);
        double end_time = GetWallTime();
        cout << "[Iteration "<< i+1  << " Time :" << (end_time-start_time) << endl;
        
        delete matrix_test;
        //delete matrix_multiply;
    }
    */
    ComplexMatrix* matrix_test = new ComplexMatrix(300,300);
    matrix_test->RandomMatrix();
    matrix_multiply = matrix_test->MultiplyToMatrix(matrix_test);
    matrix_multiply->PrintMatrix();
    ComplexMatrix* matrix_multiply_new = new ComplexMatrix(300,300);
    Complex tmp_value = 0.0;
    for(int i=0;i<matrix_test->get_row();++i)
        for(int k=0;k<matrix_test->get_column();++k)
        {
            tmp_value = 0.0;
            for(int j=0;j<matrix_test->get_column();++j)
            {
                tmp_value += matrix_test->get_matrix_element(i,j) * matrix_test->get_matrix_element(j,k);
            }
            matrix_multiply_new->set_matrix_element(i, k, tmp_value);
        }
    
    cout << "=================================" << endl;
    cout << "=================================" << endl;
    cout << "=================================" << endl;
    cout << "=================================" << endl;
    cout << "=================================" << endl;
    matrix_multiply_new->PrintMatrix();
    //for(int i=0;i<matrix_multiply->get_row();++i)
    //    for(int j=0;j<matrix_multiply->get_column();++j)
    

    delete matrix_test;
    
    cout << endl;
    cout << "7. ChangeMatrix TransposeMatrix ReshapeMatrix ExpanMatrix" << endl;
    matrix->PrintMatrix();
    
    cout << "Transpose matrix: " << endl;
    ComplexMatrix* matrix_transpose = matrix->TransposeMatrix();
    matrix_transpose->PrintMatrix();

    cout << "Reshape matrix: (2,2)->(1,4)" << endl;
    ComplexMatrix* matrix_reshape = matrix->ReshapeMatrix(1,4);
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
    matrix = new ComplexMatrix(2,3);
    for(int i=0;i<2;++i) for(int j=0;j<3;++j)
        matrix->set_matrix_element(i,j,i*3+j+1);
    matrix->PrintMatrix();
    ComplexMatrix *left_matrix, *right_matrix;
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

    cout << "SVD matrix" << endl;
    matrix_test = new ComplexMatrix(200,100);
    matrix_test->RandomMatrix();
    ComplexMatrix* matrix_test_old = new ComplexMatrix(matrix_test);
    //matrix_test->PrintMatrix();
    matrix_test->SVDMatrix(left_matrix, right_matrix, sigular_value, sigular_dim);
    ComplexMatrix* S = new ComplexMatrix(left_matrix->get_column(), right_matrix->get_row());
    for(int i=0;i<min(left_matrix->get_column(), right_matrix->get_row());++i)
    {
        S->set_matrix_element(i,i,sigular_value[i]);
    }
    ComplexMatrix* matrix_test_new;
    matrix_test_new = left_matrix->MultiplyToMatrix(S)->MultiplyToMatrix(right_matrix);
    //matrix_test_new->PrintMatrix();
    //
    for(int i=0;i<matrix_test->get_row();++i)
        for(int j=0;j<matrix_test->get_column();++j)
        {
            if(Norm(matrix_test_old->get_matrix_element(i,j) - matrix_test_new->get_matrix_element(i,j)) > 1E-8)
                error("error svd");
                //cout << Norm(matrix_test->get_matrix_element(i,j) - matrix_test_new->get_matrix_element(i,j)) << endl;
        }

    delete left_matrix, right_matrix;
    delete[] sigular_value;
    delete matrix;
    
    /* 
    cout << "diag matrix" << endl;
    double *eigenvalue=new double[30];
    int vector_dim=30;
    matrix_test = new ComplexMatrix(30,30);
    matrix_test->RandomMatrix();
    ComplexSymMatrixDiag(matrix_test->get_matrix_element(), eigenvalue, vector_dim);
    matrix_test->PrintMatrix();
    */

        
    cout << endl;
    cout << "9. ParallelMatrix" << endl;
    cout << "Left Parallel" << endl;
    matrix = new ComplexMatrix(2,3);
    matrix->set_matrix_element(0,0,1);
    matrix->set_matrix_element(0,1,3);
    matrix->set_matrix_element(0,2,1);
    matrix->set_matrix_element(1,0,2);
    matrix->set_matrix_element(1,1,4);
    matrix->set_matrix_element(1,2,2);
    matrix->PrintMatrix();
    
    int num_unparallel;
    int* position_unparallel;
    ComplexMatrix* transfer_tensor;
    matrix->ParallelMatrix(0, num_unparallel, position_unparallel, 
            transfer_tensor);
    cout << "num_unparallel: " << num_unparallel << endl;
    cout << "position_unparallel: "  << endl;
    for(int i=0;i<num_unparallel;++i)
    {
      cout << "[" << i << "] " << position_unparallel[i] << endl;
    }
    cout << "transfer_tensor: " << endl;
    transfer_tensor->PrintMatrix();

    cout << "Right Parallel" << endl;
    matrix = new ComplexMatrix(3,2);
    matrix->set_matrix_element(0,0,1);
    matrix->set_matrix_element(0,1,2);
    matrix->set_matrix_element(1,0,3);
    matrix->set_matrix_element(1,1,4);
    matrix->set_matrix_element(2,0,1);
    matrix->set_matrix_element(2,1,2);
    matrix->PrintMatrix();
    matrix->ParallelMatrix(1, num_unparallel, position_unparallel, 
            transfer_tensor);
    cout << "num_unparallel: " << num_unparallel << endl;
    cout << "position_unparallel: "  << endl;
    for(int i=0;i<num_unparallel;++i)
    {
      cout << "[" << i << "] " << position_unparallel[i] << endl;
    }
    cout << "transfer_tensor: " << endl;
    transfer_tensor->PrintMatrix();
    delete matrix;
    
    cout << endl;
    cout << "10. MatrixKronProduct" << endl;
    matrix = new ComplexMatrix(2,2);
    matrix->set_matrix_element(0,0,1);
    matrix->set_matrix_element(0,1,2);
    matrix->set_matrix_element(1,0,3);
    matrix->set_matrix_element(1,1,4);
    matrix->PrintMatrix();
    matrix->MatrixKronProduct(matrix)->PrintMatrix();

    delete matrix;
    
    return 0;
}
