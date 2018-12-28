#include "util/general.h"
#include "util/general.ih"

void ReorderRelevantArray(int num_quantum, int num_table, 
        int** quantum_table, int* index)
{
    int **tmp_quantum_table;
    tmp_quantum_table = new int* [num_table];
    for(int i=0;i<num_table;++i)
        tmp_quantum_table[i] = new int[num_quantum];
    // reorder
    for(int i=0;i<num_table;++i)
        for(int j=0;j<num_quantum;++j)
            tmp_quantum_table[index[i]][j] = quantum_table[i][j];
    // quantum_table = tmp_quantum_table
    for(int i=0;i<num_table;++i)
        for(int j=0;j<num_quantum;++j)
            quantum_table[i][j] = tmp_quantum_table[i][j];

    for(int i=0;i<num_table;++i) delete[] tmp_quantum_table[i];
    delete[] tmp_quantum_table;
}
