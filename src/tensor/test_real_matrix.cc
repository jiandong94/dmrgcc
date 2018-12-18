#include "tensor/real_matrix.h"

//class RealMatrix;
int main()
{
    // constructor
    RealMatrix* matrix1 = new RealMatrix();
    delete matrix1;

    // constructor
    int m=2,n=2;
    cout << "matrix2:" << endl;
    RealMatrix* matrix2 = new RealMatrix(m,m);
    for(int i=0;i<m;++i) 
        for(int j=0;j<m;++j)
            matrix2->set_matrix_element(i,j,i*m+j+1);

    // PrintMatrix
    matrix2->PrintMatrix();

    //get
    int row = matrix2->get_row();
    int column = matrix2->get_column();
    int total_element_num = matrix2->get_total_element_num();
    int position_0_0 = matrix2->get_matrix_element(0,0);
    double* matrix_element = matrix2->get_matrix_element();
    cout << "row = " << row << endl;
    cout << "column = " << column << endl;
    cout << "total_element_num = " << total_element_num << endl;
    cout << "position_0_0 =  " << position_0_0 << endl;
    cout << "position_0_0 =  " << matrix_element[0] << endl;


    cout << "============================================================" << endl;
    // RandomMatrix
    // ClearMatrix
    RealMatrix* matrix3 = new RealMatrix(m,n);
    cout << "(matrix3)Random Matrix:" << endl;
    matrix3->RandomMatrix();
    matrix3->PrintMatrix();

    RealMatrix* matrix4 = new RealMatrix(m,n);
    cout << "(matrix4)Clear Matrix:" << endl;
    matrix4->ClearMatrix();
    matrix4->PrintMatrix();

    //matrix3->ResetMatrix();
    //cout << "after ResetMatrix, total_element_num = " <<matrix3->get_total_element_num() << endl;


    // add element
    // add matrix
    cout << "matrix4(0,0)+1 AddToMatrixElement:" << endl;
    matrix4->AddToMatrixElement(0,0,1);
    matrix4->PrintMatrix();

    cout << "matrix4+matrix3 AddToMatrix:" << endl;
    matrix4->AddToMatrix(matrix3);
    matrix4->PrintMatrix();

    cout << "matrix4*2 MultiplyToScalar:" << endl;
    matrix4->MultiplyToScalar(2);
    matrix4->PrintMatrix();

    cout << "matrix4*matrix3 MatrixElementProduct:" << endl;
    matrix4->MatrixElementProduct(matrix3);
    matrix4->PrintMatrix();

    cout << "matrix2*matrix2 MultiplyToMatrix:" << endl;
    RealMatrix* matrix5;
    matrix5 = new RealMatrix(matrix2->MultiplyToMatrix(matrix2));
    matrix5->PrintMatrix();
    
    delete matrix3;
    delete matrix4;
    delete matrix5;
    // test speed
    cout << "============================================================" << endl;
    cout << "Test MKL speed" << endl;
    for(int i=0;i<10;++i)
    {
        RealMatrix* matrix6 = new RealMatrix(3000,3000);
        matrix6->RandomMatrix();
        double start_time = GetWallTime();
        RealMatrix* matrix7 = matrix6->MultiplyToMatrix(matrix6);
        double end_time = GetWallTime();
        cout << "[Iteration "<< i+1  << " Time :" << (end_time-start_time) << endl;
        
        delete matrix6;
        delete matrix7;
    }
    
    // transpose matrix
    // reshape matrix
    cout << "" << endl;
    cout << "test transpose matrix" << endl;
    RealMatrix* transpose_matrix = matrix2->TransposeMatrix();
    transpose_matrix->PrintMatrix();


    cout << "test reshape matrix" << endl;
    RealMatrix* reshape_matrix = matrix2->ReshapeMatrix(1,4);
    reshape_matrix->PrintMatrix();

    delete transpose_matrix, reshape_matrix;



    // combine matrix
    // svdmatrix
    cout << "" << endl;
    cout << "test combine matrix" << endl;
    matrix1 = new RealMatrix(matrix2);
    matrix1->PrintMatrix();
    matrix2->PrintMatrix();

    matrix1->CombineMatrix(1, matrix2);
    matrix1->PrintMatrix();


    cout << "test svd" << endl;
    RealMatrix *left_matrix, *right_matrix;
    double* sigular_value;
    int sigular_dim;
    matrix2->SVDMatrix(left_matrix, right_matrix, sigular_value, sigular_dim);
    left_matrix->PrintMatrix();
    right_matrix->PrintMatrix();
    cout << "sigular dim: " << sigular_dim << endl;
    cout << "sigular value: ";
    for(int i=0;i<sigular_dim;++i) 
        cout << sigular_value[i] << " ";
    cout << endl;

    
    delete left_matrix, right_matrix;
    delete[] sigular_value;
    delete matrix1, matrix2;

    // write read
    char* test_matrix = "test_matrix.dat";
    matrix1 = new RealMatrix(5,5);
    matrix1->RandomMatrix();
    matrix1->PrintMatrix();

    matrix1->WriteMatrix(test_matrix);
    matrix2 = new RealMatrix();
    matrix2->ReadMatrix(test_matrix);
    matrix2->PrintMatrix();

    delete matrix1, matrix2;


    return 0;
}
