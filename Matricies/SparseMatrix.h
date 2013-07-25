//
//  SparseSparseMatrix.h
//  Graph Matching
//
//  Created by Ali Hajimirza on 7/9/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _SparseMatrix_h
#define _SparseMatrix_h

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include "SparseElement.h"
#include "MatrixExceptions.h"

#ifdef __linux__
#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"
#endif

template <typename T>
class SparseMatrix;

template <typename T>
std::ostream& operator<< (std::ostream&, const SparseMatrix<T>&);

template <typename T>
class SparseMatrix
{
private:
    const static int _DEFAULT_MATRIX_SIZE;
    const static T _DEFAULT_MATRIX_ENTRY;
    size_t _getArrSize() const;
    void _copy(const SparseMatrix<T>&);
    bool _hasEdge(int,int);
    
protected:
    int _rows;
    int _cols;
    std::unordered_map <int, T > _edges;   
public:
    /**************
     *Constructors*
     **************/
    SparseMatrix();
    SparseMatrix(int, int);
    SparseMatrix(const std::string&);
    SparseMatrix(const SparseMatrix<T>&);
    
    /************
     *Destructor*
     ************/
    virtual ~SparseMatrix();
    
    /***********
     *ACCESSORS*
     ***********/
    bool isSquare() const;
    bool isSymmetric() const;
    int getNumberOfRows();
    int getSparseFormSize();
    int getNumberOfColumns();
    std::vector<T> getSumOfRows();
    std::vector<T> getTopEigenVector();
    std::vector<int> getNeighbors(int vertex);
    //sparse_SparseMatrix_element<T>** getSparseForm();
    SparseMatrix<T> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B);
    
    /**********
     *MUTATORS*
     **********/
    void insert(int, int, T);

    /**********
     *OPERATORS*
     **********/
    T operator()(int,int) const;
    SparseMatrix<T> kron(const SparseMatrix<T>&); /* WORKS SHOULD BE CHANGED */
    bool operator==(const SparseMatrix<T>&);
    SparseMatrix<T>& operator- (SparseMatrix<T>& other_matrix);
    SparseMatrix<T>& operator= (const SparseMatrix<T>&);
    SparseMatrix<T> diagonalVectorTimesSparseMatrix(const std::vector<T>&);
    SparseMatrix<T> SparseMatrixTimesDiagonalVector(const std::vector<T>&);
    friend std::ostream& operator<< <> (std::ostream& stream, const SparseMatrix<T>& SparseMatrix);

    /* WILL IMPLEMENT IF I HAD TIME
     SparseMatrix<T> operator* (const SparseMatrix<T>&);
     SparseMatrix<T> operator+ (const SparseMatrix<T>&);
     SparseMatrix<T> operator- (const SparseMatrix<T>&);
     SparseMatrix<T> operator* (T);
     SparseMatrix<T> operator+ (T);
     SparseMatrix<T> operator- (T);
     void operator*= (const SparseMatrix<T>&);
     void operator+= (const SparseMatrix<T>&);
     void operator-= (const SparseMatrix<T>&);
     void operator*= (T);
     void operator+= (T);
     void operator-= (T);
     */
    
};
//==========================================================CONSTANTS============================================================
template <typename T>
const int SparseMatrix<T>::_DEFAULT_MATRIX_SIZE = 1;
template <typename T>
const T SparseMatrix<T>::_DEFAULT_MATRIX_ENTRY = 1;
//==========================================================CONSTRUCTORS============================================================

template <typename T>
inline SparseMatrix<T>::SparseMatrix()
{
    this->_rows = _DEFAULT_MATRIX_SIZE;
    this->_cols = _DEFAULT_MATRIX_SIZE;
}

template<typename T>
inline SparseMatrix<T>::SparseMatrix(const std::string& file_path)
{
    int tmp_i;
    int tmp_j;
    std::ifstream file_reader;
    file_reader.open(file_path.c_str());
    
    if(file_reader.fail())
    {
        file_reader.close();
        throw FileDoesNotExistException(file_path.c_str());
    }
    
    file_reader >> this->_rows;
    file_reader >> this->_cols;
    file_reader >> tmp_i;       //skip the number of lines entry
    
    while (!file_reader.eof())
    {
        file_reader >> tmp_i;
        file_reader >> tmp_j;
        this->insert(tmp_i -1, tmp_j -1 , 1);
    }
    file_reader.close();
}

template <typename T>
inline SparseMatrix<T>::SparseMatrix(int rows, int cols)
{
    this->_rows = rows;
    this->_cols = cols;
}

template <typename T>
inline SparseMatrix<T>::SparseMatrix(const SparseMatrix<T>& SparseMatrix)
{
    _copy(SparseMatrix);
}

//==========================================================DESTRUCTOR==============================================================
template <typename T>
inline SparseMatrix<T>::~SparseMatrix()
{
}

//===========================================================ACCESSORS===============================================================
template <typename T>
inline int SparseMatrix<T>::getNumberOfRows()
{
    return this->_rows;
}

template <typename T>
inline int SparseMatrix<T>::getNumberOfColumns()
{
    return this->_cols;
}

template <typename T>
inline bool SparseMatrix<T>::isSquare() const
{
    return (this->_rows == this->_cols);
}

template <typename T>
inline bool SparseMatrix<T>::isSymmetric() const ////////////////////////////TODO
{
    // //sym. SparseMatrix has to be square
    // if (this->_rows != this->_cols)
    // {
    //     return false;
    // }
    
    // //chaking for entries to be equal
    // for(int i = 0; i < this->_rows; i++)
    // {
    //     for(int j = 0; j < this->_cols; j++)
    //     {
    //         if (this->_edges[(i * this->_cols) + j] != this->_edges[(j * this->_cols) + i])
    //         {
    //             return false;
    //         }
    //     }
    // }
    
    // return true;
}

//template <typename T>
//sparse_SparseMatrix_element<T>** SparseMatrix<T>::getSparseForm(){
//    
//    if(sparse_form != NULL)
//    {
//        return this->sparse_form;
//    }
//    
//    for(int i = 0; i < this->_getArrSize(); i++)
//    {
//        if(this->_edges[i] != 0)
//        {
//            sparse_form_size++;
//        }
//    }
//    
//    sparse_form = new sparse_SparseMatrix_element<T>* [sparse_form_size];
//    int counter=0;
//    
//    for(int i = 0; i < this->_getArrSize(); i++)
//    {
//        if(this->_edges[i] != 0)
//        {
//            sparse_SparseMatrix_element<T> *ele= new sparse_SparseMatrix_element<T>;
//            ele->row_index = i / this->_rows;
//            ele->col_index = i % this->_rows;
//            ele->value = this->_edges[i];
//            sparse_form[counter]=ele;
//            counter++;
//        }
//    }
//    
//    return this->sparse_form;
//}

template <typename T>
inline SparseMatrix<T> SparseMatrix<T>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B)
{
    int num_in_A = 0;
    for (int i=0; i< vec_A.size(); i++)
    {
        if (vec_A[i] == 1)
        {
            num_in_A++;
        }
    }
    int num_in_B = 0;
    for (int i=0; i < vec_B.size(); i++)
    {
        if( vec_B[i] == 1)
        {
            num_in_B++;
        }
    }
    //Initializing and allocating the product SparseMatrix
    SparseMatrix<T> res_SparseMatrix(num_in_A, num_in_B);
    
    int counter = 0;
    
    for (int i=0; i < vec_A.size(); i++)
    {
        for(int j=0; j < vec_B.size(); j++)
        {
            if ( vec_A[i] == 1 && vec_B[j] ==1)
            {
                res_SparseMatrix._edges[counter] = (*this)(i,j);
                counter++;
            }
        }
    }
    return res_SparseMatrix;
}



template <typename T>
inline std::vector<int> SparseMatrix<T>::getNeighbors(int vertex)
{
    std::vector<int> neighbors;
    
    for(int i = 0; i < this->getNumberOfRows(); i++)
    {
        if(this->_edges[(i * this->_cols) + vertex] == 1)
        {
            neighbors.push_back(i);
        }
    }
    
    return neighbors;
}

template <typename T>
inline std::vector<T> SparseMatrix<T>::getTopEigenVector(){
    
    //    int arr_size= (0.5 * this->_rows * (this->_rows+1));
    //    double nzval[arr_size];
    //    int counter = 0;
    //
    //    for(int i = 0; i < this->_rows; i++)
    //    {
    //        for(int j = i; j < this->_rows; j++)
    //        {
    //            nzval[counter] = (*this)(i,j);
    //            counter++;
    //        }
    //    }
    //
    //    ARdsSymSparseMatrix<T> ARSparseMatrix(this->_rows,nzval,'L');
    //    ARluSymStdEig<T> eigProb(1, ARSparseMatrix, "LM", 10);
    //    eigProb.FindEigenvectors();
    //    std::vector<T>* eigenVec = new std::vector<T>(eigProb.GetN());
    //
    //    for (int i=0; i < eigProb.GetN() ; i++)
    //    {
    //        eigenVec[i] = eigProb.Eigenvector(0,i);
    //    }
    //
    //    return eigenVec;
}


/*
 * returns a vector of type T where each entry corresponds to sum of entries in a SparseMatrix row.
 * Delete this object after usage.
 * @return std::vector<T>
 */
template <typename T>
inline std::vector<T> SparseMatrix<T>::getSumOfRows()
{
    std::vector<T>* sum_vector= new std::vector<T>(this->_rows);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        (*sum_vector)[(i / this->_rows)] += this->_edges[i];
    }
    
    return sum_vector;
}

//===========================================================MUTATORS================================================================
template <typename T>
inline void SparseMatrix<T>::insert(int i, int j, T value)
{
    //if the edge exist =>delete the node
    if (_hasEdge(i,j))
    {
            this->_edges.erase((i*this->_cols) + j);
    }
    //if the evalue is not 0 =>insert it
    if (value != 0)
    {
        this->_edges[(i*this->_cols) + j] = value;
    }
}


//==========================================================OPERATORS================================================================
template <typename T>
inline SparseMatrix<T> SparseMatrix<T>::kron(const SparseMatrix<T>& matrix)
{
    // checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
       throw NotASquareMatrixException();
    }
    
   
    int prod_size = this->_rows * matrix._rows;
    SparseMatrix<T> prod_matrix(prod_size, prod_size);

    int mat1_i;
    int mat1_j;
    int mat2_i;
    int mat2_j;
     for(auto i_it = this->_edges.begin(); i_it != this->_edges.end(); ++i_it)
    {
        for(auto j_it = matrix._edges.begin(); j_it != matrix._edges.end(); ++j_it)
        {
            mat1_i = i_it->first/this->_cols;
            mat1_j = i_it->first%this->_cols;
            mat2_i = j_it->first/matrix._cols;
            mat2_j = j_it->first%matrix._cols;

            prod_matrix._edges[(mat1_i * matrix._cols + mat2_i)*prod_size+ (mat1_j * matrix._cols + mat2_j)] = (i_it->second) * (j_it->second);
        }
    }
    
    return prod_matrix;
}

template <typename T>
inline SparseMatrix<T> SparseMatrix<T>::diagonalVectorTimesSparseMatrix(const std::vector<T>& vec)
{
    if(this->_rows != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    SparseMatrix<T> ret_SparseMatrix (*this);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        ret_SparseMatrix._edges[i] = vec[i/this->_rows] * this->_edges[i];
    }
    
    return ret_SparseMatrix;
}

template <typename T>
inline SparseMatrix<T> SparseMatrix<T>::SparseMatrixTimesDiagonalVector(const std::vector<T>& vec)
{
    if(this->_cols != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    SparseMatrix<T> ret_SparseMatrix(*this);
    for(int i = 0; i < this->_getArrSize(); )
    {
        for(int j = 0; j < vec.size(); j++)
        {
            ret_SparseMatrix._edges[i] = this->_edges[i] * vec[j];
            i++;
        }
    }
    
    return ret_SparseMatrix;
}
/////////////////////////////////this point
template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const SparseMatrix<T>& matrix)
{
    stream<< "Size: " << matrix._rows << "*" << matrix._cols << '\n';
    for (int i = 0; i < matrix._rows; i++)
    {
        for (int j = 0; j < matrix._cols; j++)
        {
            stream << matrix(i, j) << ' ';
        }
        stream << "\n";
    }
    stream << "\n\n\n";
    return stream;
}

template <typename T>
inline T SparseMatrix<T>::operator()(int i, int j) const
{
    if ( (i >= _rows) || (i < 0) || (j >= _cols) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    typename std::unordered_map<int,T>::const_iterator find_res = this->_edges.find((i * this->_cols) + j);
    if ( find_res == this->_edges.end() )
        return 0;
    else
        return find_res->second;
}


template <typename T>
inline SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix<T>& rhs)
{
    _copy(rhs);
}

template <typename T>
inline bool SparseMatrix<T>::operator==(const SparseMatrix<T>& rhs)
{
    // checking for dimension equality
    if ((this->_cols != rhs._cols) || (this->_rows != rhs._rows) || (this->_getArrSize() != rhs._getArrSize()))
    {
        return false;
    }
    
    //checking for entry equality
    for (int i = 0; i < this->_getArrSize(); i++)
    {
        for( int j = 0; j< this->_cols; j++)
        {
            if ((*this)(i,j) != rhs(i,j))
                return false;
        }
    }
    
    return true;
}

//===========================================================PRIVATE=================================================================
template <typename T>
inline void SparseMatrix<T>::_copy(const SparseMatrix<T>& rhs)
{
    this->_rows = rhs._rows;
    this->_cols = rhs._cols;
    this->_edges = std::unordered_map<int, T>(rhs._edges);
}

template<typename T>
inline size_t SparseMatrix<T>::_getArrSize() const
{
    return this->_edges.size();
}

template <typename T>
inline bool SparseMatrix<T>::_hasEdge(int i, int j)
{
    return (this->_edges.find((i * this->_cols) + j) != this->_edges.end());
}


//===================================================================================================================================
    

#endif
