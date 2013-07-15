//
//  SymMatrix.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/24/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _SymMatrix_h
#define _SymMatrix_h


#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <cstring>
#include <limits>
#include "Matrix.h"
#include "MatrixExceptions.h"

#ifdef __linux__
#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"
#endif

template <typename DT>
class SymMatrix;

template <typename DT>
std::ostream& operator<< (std::ostream&, SymMatrix<DT>&);

template <typename DT>
class SymMatrix
{
private:
    const static int _DEFAULT_MATRIX_SIZE = 1;
    size_t _getArrSize() const;
    void _initilalizeMatrix();
    void _copy(const SymMatrix<DT>&);
    
protected:
    int _size;
    DT* _edges;

public:
    
    /**************
     *Constructors*
     **************/
    SymMatrix();
    SymMatrix(int size);
    SymMatrix(const SymMatrix<DT>&);
    SymMatrix(const std::string&);
    
    /************
     *Destructor*
     ************/
    virtual ~SymMatrix();
    
    /***********
     *ACCESSORS*
     ***********/
    int getSize();
    int getSparseFormSize();
    std::vector<DT> getSumOfRows();
    std::vector<DT> getTopEigenVector();
    std::vector<int> getNeighbors(int vertex);
    Matrix<DT> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B);
    
    /**********
     *MUTATORS*
     **********/
    void insert(int, int, DT value);
    
    /**********
     *OPERATORS*
     **********/
    DT operator()(int,int);
    void operator()(int,int,DT);
    bool operator==(const SymMatrix<DT>&);
    SymMatrix<DT> kron(const SymMatrix<DT>&); /* WORKS SHOULD BE CHANGED */
    SymMatrix<DT>& operator= (const SymMatrix<DT>&);
    SymMatrix<DT> diagonalVectorTimesMatrix(const std::vector<DT>&);
    SymMatrix<DT> matrixTimesDiagonalVector(const std::vector<DT>&);
    friend std::ostream& operator<< <> (std::ostream& stream, const Matrix<DT>& matrix);
    
    
    /* WILL IMPLEMENT IF I HAD TIME AND IF THE OPERATION MADE SENSE FOR SYMMATRIX
     SymMatrix<DT> operator* (const SymMatrix<DT>&);
     SymMatrix<DT> operator+ (const SymMatrix<DT>&);
     SymMatrix<DT> operator- (const SymMatrix<DT>&);
     SymMatrix<DT> operator* (DT);
     SymMatrix<DT> operator+ (DT);
     SymMatrix<DT> operator- (DT);
     void operator*= (const SymMatrix<DT>&);
     void operator+= (const SymMatrix<DT>&);
     void operator-= (const SymMatrix<DT>&);
     void operator*= (DT);
     void operator+= (DT);
     void operator-= (DT);
     */
};

//==========================================================CONSTRUCTORS============================================================

template <typename DT>
inline SymMatrix<DT>::SymMatrix()
{
    this->_size = _DEFAULT_MATRIX_SIZE;
    _initilalizeMatrix();
}

template<typename DT>
inline SymMatrix<DT>::SymMatrix(const std::string& file_path)
{
    int tmp_x;
    int tmp_y;
    std::ifstream file_reader;
    file_reader.open(file_path.c_str());
    
    if(file_reader.fail())
    {
        file_reader.close();
        throw FileDoesNotExistException(file_path.c_str());
    }
    
    file_reader >> this->_size;
    file_reader >> tmp_x;
    
    if (this->_size != tmp_x)
    {
        file_reader.close();
        throw NotASquareMatrixException();
    }
    
    _initilalizeMatrix();
    file_reader >> tmp_x;      //skip the line number
    
    while (!file_reader.eof())
    {
        file_reader >> tmp_x;
        file_reader >> tmp_y;
        this->insert(tmp_x - 1 ,tmp_y - 1, 1);
    }
    
    file_reader.close();
}

template <typename DT>
inline SymMatrix<DT>::SymMatrix(int size)
{
    this->_size = size;
    _initilalizeMatrix();
}

template <typename DT>
inline SymMatrix<DT>::SymMatrix(const SymMatrix<DT>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
template <typename DT>
inline SymMatrix<DT>::~SymMatrix()
{
    if (_edges != NULL)
    {
        delete [] _edges;
    }
}

//===========================================================ACCESSORS===============================================================
template <typename DT>
inline int SymMatrix<DT>::getSize()
{
    return _size;
}

template <typename DT> ///////Needs to be changed
inline SymMatrix<DT> SymMatrix<DT>::kron(const SymMatrix<DT>& matrix)
{
    // also checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    //Initializing and allocating the product matrix
    int prod_size = this->_size * matrix._size;
    SymMatrix<DT> prod_matrix(prod_size);
    
    /*
     *  Calculating the kronecker product:
     *  The indices of the product matrix is calculated by:
     *      i = (i_outer*size) + i_inner
     *      j = (j_outer*size) + j_inner
     */
    for (int i_outer = 0; i_outer < this->_size; i_outer++)
    {
        for (int j_outer=0; j_outer < this->_size; j_outer++)
        {
            for(int i_inner=0; i_inner < matrix._size; i_inner++)
            {
                for(int j_inner=0; j_inner < matrix._size; j_inner++)
                {
                    prod_matrix.insert((i_outer*matrix._size) + i_inner,
                                        (j_outer*matrix._size)+j_inner,
                                        (*this)(i_outer, j_outer) * matrix(i_inner, j_inner));
                }
            }
        }
    }
    
    return prod_matrix;
}

template <typename DT> //NEEDS REVIEW
inline Matrix<DT> SymMatrix<DT>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B)
{
    int num_in_A = 0;
    for (int i=0; i< vec_A.size(); i++)
    {
        if (vec_A[i] != 0)
        {
            num_in_A++;
        }
    }
    int num_in_B = 0;
    for (int i=0; i < vec_B.size(); i++)
    {
        if( vec_B[i] != 0)
        {
            num_in_B++;
        }
    }
    //Initializing and allocating the product matrix
    Matrix<DT> res_matrix(num_in_A, num_in_B);
    
    int counter = 0;
    
    for (int i=0; i< vec_A.size(); i++)
    {
        for(int j=0; j< vec_B.size(); j++)
        {
            if ( vec_A[i] == 1 && vec_B[j] ==1)
            {
                res_matrix.insert(counter/num_in_B, counter%num_in_B, (*this)(i, j));
                counter++;
            }
        }
    }
    return *res_matrix;
}



template <typename DT> //NEEDS REVIEW
inline std::vector<int> SymMatrix<DT>::getNeighbors(int vertex)
{
    std::vector<int> neighbors;
    
    for(int i=0; i < this->size(); i++)
    {
        if( (*this)(i,vertex) == 1 )
        {
            neighbors.push_back(i);
        }
    }
    
    return neighbors;
}

template <typename DT>
inline std::vector<DT> SymMatrix<DT>::getTopEigenVector()
{
#ifdef __linux__
    ARdsSymMatrix<DT> ARMatrix(this->_size, this->_edges, 'L');
    ARluSymStdEig<DT> eigProb(1, ARMatrix, "LM", 10);
    eigProb.FindEigenvectors();
    std::vector<DT> eigen_vec = new std::vector<DT>(eigProb.GetN());
     
    for (int i=0; i < eigProb.GetN() ; i++)
    {
        eigen_vec[i] = eigProb.Eigenvector(0,i);
    }
    
    return eigen_vec;
#endif
}


//===========================================================MUTATORS================================================================

/*
 * returns an array of floats where each array element corresponds to sum of entries in a matrix row
 * @pram: pointer to 2-d array of floats
 * @pram: number of rows in graph
 * @pram: number of columns in graph
 */
template <typename DT>
inline std::vector<DT> SymMatrix<DT>::getSumOfRows()
{
    std::vector<DT> result_vec(this->_size);
    
    for (int i = 0 ; i < this->_size; i++)
    {
        for (int j = 0; j < this->_size; j++)
        {
            (*result_vec)[i] = (*result_vec)[i] + (*this)(i, j);
        }
    }
    
    return result_vec;
}

template <typename DT>
inline void SymMatrix<DT>::insert(int i, int j, DT value)
{
    if ( (i >= _size) || (i < 0) || (j >= _size) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    
    if(i<=j)
    {
        this->_edges[i+size_t(j)*(j+1)/2] = value ;
    }
    else
    {
        this->_edges[j+size_t(i)*(i+1)/2] = value;
    }
}



//==========================================================OPERATIONS================================================================
template <typename DT>
inline SymMatrix<DT> SymMatrix<DT>::diagonalVectorTimesMatrix(const std::vector<DT>& vec)
{
    if(this->_rows != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    Matrix<DT> ret_matrix(*this);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        ret_matrix._edges[i] = vec[i/this->_rows] * this->_edges[i];
    }
    
    return ret_matrix;
}

template <typename DT>
inline SymMatrix<DT> SymMatrix<DT>::matrixTimesDiagonalVector(const std::vector<DT>& vec)
{
    if(this->_cols != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    Matrix<DT> ret_matrix(*this);
    for(int i = 0; i < this->_getArrSize(); )
    {
        for(int j = 0; j < vec.size(); j++)
        {
            ret_matrix._edges[i] = this->_edges[i] * vec[j];
            i++;
        }
    }
    
    return ret_matrix;
}

//==========================================================OPERATORS================================================================
template <typename DT>
inline std::ostream& operator<<(std::ostream& stream, SymMatrix<DT>& matrix)
{
    stream<< "Size: " << matrix._size << "*" << matrix._size << '\n';
    for (int i = 0; i < matrix._size; i++)
    {
        for( int j = 0; j< matrix._size; j++)
        {
            stream << matrix(i, j) << ' ';
        }
        stream << "\n";
    }
    stream << "\n\n\n";
    return stream;
}

template <typename DT>
inline DT SymMatrix<DT>::operator()(int i, int j)
{
    if ( (i >= _size) || (i < 0) || (j >= _size) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    
    if(i<=j)
    {
        return this->_edges[i+size_t(j)*(j+1)/2];
    }
    else
    {
        return this->_edges[j+size_t(i)*(i+1)/2];
    }
}
    
template <typename DT>
inline void SymMatrix<DT>::operator()(int i, int j, DT value)
{
    if ( (i >= _size) || (i < 0) || (j >= _size) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    
    if(i<=j)
    {
        this->_edges[i+size_t(j)*(j+1)/2] = value;
    }
    else
    {
        this->_edges[j+size_t(i)*(i+1)/2] = value;
    }
}

template <typename DT>
inline SymMatrix<DT>& SymMatrix<DT>::operator=(const SymMatrix<DT>& matrix)
{
    delete [] this->_edges;
    _copy(matrix);
}

//===========================================================PRIVATE=================================================================
template <typename DT>
inline void SymMatrix<DT>::_copy(const SymMatrix<DT>& matrix)
{
    this->_size = matrix._size;
    !_initilalizeMatrix();
  
    for(int i=0; i < this->_getArrSize() ; i++)
    {
        this->_edges[i] = matrix._edges[i];
    }
}

template <typename DT>
inline void SymMatrix<DT>::_initilalizeMatrix()
{
    this->sparse_form_size = 0;
    this->sparse_form = NULL;
    this->_edges = new DT[this->_getArrSize()];
    
    if (this->_edges == NULL)           //checking for memory issues.
    {
        throw OutOfMemoryException();
    }
    
    for ( int i = 0; i < this->_getArrSize(); i++)
    {
        this->_edges[i] = 0;
    }
}
    
template<typename DT>
inline size_t SymMatrix<DT>::_getArrSize() const
{
    return (int)(0.5 * this->_size * (this->_size + 1));
}

//===================================================================================================================================

#endif
