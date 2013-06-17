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
#include "util.h"
#include <vector>

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
    std::cout << *kron_prod << std::endl;
    
     graph_con_com(kron_prod, kron_prod->getNumberOfColumns(),&vertex_stack);
//    
////    for(int i=0; i < kron_prod->getNumberOfColumns(); i++ )
////    {
//        int* comp_mask=component_mask(vertices, 0, kron_prod->getNumberOfColumns());
//        
//        // Temp part should be replaced
//        std::vector<int> temp_vec(kron_prod->getNumberOfRows());
//        for(int j=0; j < kron_prod->getNumberOfRows() ;j++ )
//        {
//            temp_vec[j] = comp_mask[j];
//            std::cout << temp_vec[j] << ", ";
//        }
//        //temp ends
//        
//        SparseMatrix<DT>* L = kron_prod->getScatteredSelection(temp_vec,temp_vec);
//        std::cout << *L << std::endl;
//        DT* sum = L->sum_rows();
////        std::vector <DT> D_neg1 (kron_prod->getNumber);
////        std::vector <DT> D_0pt5 (kron_prod->getNumber);
////        std::vector <DT> D_neg0pt5 (kron_prod->getNumber);
//        
//        DT* D_neg1  = new DT[kron_prod->getNumberOfRows()];
//        DT* D_0pt5  = new DT[kron_prod->getNumberOfRows()];
//        DT* D_neg0pt5  = new DT[kron_prod->getNumberOfRows()];
//        
//        for(int j=0; j < kron_prod->getNumberOfRows(); j++)
//        {
//            D_neg1[j] = 1.0/sum[j];
//            D_0pt5[j] = sqrt(sum[j]);
//            D_neg0pt5[j] = 1.0/D_0pt5[j];
//        }
//        
//        //SparseMatrix<DT>* M = L->vec_times_mat(D_neg1,kron_prod->getNumberOfRows() );
//       // SparseMatrix<DT>* Ms = (L->vec_times_mat(D_0pt5,kron_prod->getNumberOfRows()))->mat_times_vec(D_0pt5,kron_prod->getNumberOfRows());
////        if(!Ms->isSymmetric())
////        {
////            throw NotASymmetricMatrixException();
////        }
//        
////    }
    
}

#endif
