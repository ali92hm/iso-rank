//
//  IsoRank.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/14/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef Sparse_Matrix_IsoRank_h
#define Sparse_Matrix_IsoRank_h

#include "SparseMatrix.h"
#include "Tarjan.h"

template <typename DT>
void isoRank(SparseMatrix<DT>& matrix_A, SparseMatrix<DT>& matrix_B)
{
    if (!matrix_A.isSquare() || !matrix_B.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    if (!matrix_A.isSymmetric() || !matrix_B.isSymmetric())
    {
        throw NotASymmetricMatrixException();
    }
    
    // Degree distribution statistics
    
    
    SparseMatrix<DT>* kron_prod = matrix_A.kron(matrix_B);
    std::stack<vertex*> vertex_stack;
    vertex* vertices= graph_con_com(kron_prod, kron_prod->getNumberOfColumns(),&vertex_stack);
    
    for(int i=0; i < kron_prod->getNumberOfColumns(); i++ )
    {
        int* comp_mask=component_mask(vertices, i, kron_prod->getNumberOfColumns());
    }
    
}

#endif
