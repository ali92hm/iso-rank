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
    
   
    vector<vertex*> * vertices = graph_con_com(kron_prod, kron_prod->getNumberOfColumns());
//    cout<< "connected comps: " <<endl;
//    for (int i=0; i< vertices->size(); i++)
//    {
//        cout << (*vertices)[i]->get_low_link() << ", ";
//    }
//    cout<< endl;
    
    for(int i=0; i < kron_prod->getNumberOfColumns(); i++ )
    {
//        cout << "Starting the for loop " << i << endl;
        vector<int>* comp_mask = component_mask(*vertices, i);
        if (comp_mask == NULL)
        {
//            cout << "comp_mask was NULL" << endl;
            continue;
        }
    
        SparseMatrix<DT>* L = kron_prod->getScatteredSelection(*comp_mask,*comp_mask);
//        cout << "Matrix L \n" << *L << endl;
        
        
        DT* sum = L->sum_rows();
        DT* D_neg1  = new DT[L->getNumberOfRows()];
        DT* D_0pt5  = new DT[L->getNumberOfRows()];
        DT* D_neg0pt5  = new DT[L->getNumberOfRows()];
        
        for(int j=0; j < L->getNumberOfRows(); j++)
        {
            D_neg1[j] = 1.0/sum[j];
            D_0pt5[j] = sqrt(sum[j]);
            D_neg0pt5[j] = 1.0/D_0pt5[j];
        }
    
        SparseMatrix<DT>* M = L->vec_times_mat(D_neg1,L->getNumberOfRows() );
        SparseMatrix<DT>* Ms = (L->vec_times_mat(D_neg0pt5,L->getNumberOfRows()))->mat_times_vec(D_neg0pt5,L->getNumberOfRows());
        if(!Ms->isSymmetric())
        {
            throw NotASymmetricMatrixException();
        }
    
        cout<< "Ms Matrix" << endl << *Ms << endl;
        double* eig = Ms->getTopEigenVector();
        

        delete comp_mask;
        delete L;
        delete [] sum;
        delete [] D_neg1;
        delete [] D_0pt5;
        delete [] D_neg0pt5;
        delete M;
        delete Ms;
    }

    typename vector<vertex*>::iterator i;
    for ( i = vertices->begin() ; i < vertices->end(); ++i )
    {
        delete * i;
    }
    delete vertices;
    delete kron_prod;
    
}

#endif
