//
//  Matrix.h
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _Matrix_h
#define _Matrix_h

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "Matrix2D.h"
#include "MatrixExceptions.h"

#ifdef __linux__
#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"
#endif

template <typename T>
class Matrix;

template <typename T>
std::ostream& operator<< (std::ostream&, const Matrix<T>&);

template <typename T>
class Matrix
{
private:
    const static int _DEFAULT_MATRIX_SIZE = 1;
    size_t _getArrSize() const;
    void _initilalizeMatrix();
    void _copy(const Matrix<T>&);
    
protected:
    int _rows;
    int _cols;
    T* _edges;
    
    /////////////////////////////////////////////
    int sparse_form_size;
    sparse_matrix_element<T>** sparse_form;
    /////////////////////////////////////////////
    
public:
    /**************
     *Constructors*
     **************/
    Matrix();
    Matrix(int, int);
    Matrix(const Matrix<T>&);
    Matrix(const std::string&);
    
    /************
     *Destructor*
     ************/
    virtual ~Matrix();
    
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
    sparse_matrix_element<T>** getSparseForm();
    Matrix<T> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B);
    
    /**********
     *MUTATORS*
     **********/
    void insert(int, int, T);
    
    /**********
     *OPERATORS*
     **********/
    T operator()(int,int) const;
    void operator()(int,int,T);
    Matrix<T> kron(const Matrix<T>&); /* WORKS SHOULD BE CHANGED */
    bool operator==(const Matrix<T>&);
    Matrix<T>& operator= (const Matrix<T>&);
    Matrix<T> diagonalVectorTimesMatrix(const std::vector<T>&);
    Matrix<T> matrixTimesDiagonalVector(const std::vector<T>&);
    friend std::ostream& operator<< <> (std::ostream& stream, const Matrix<T>& matrix);
    
    
    /* WILL IMPLEMENT IF I HAD TIME
     Matrix<T> operator* (const Matrix<T>&);
     Matrix<T> operator+ (const Matrix<T>&);
     Matrix<T> operator- (const Matrix<T>&);
     Matrix<T> operator* (T);
     Matrix<T> operator+ (T);
     Matrix<T> operator- (T);
     void operator*= (const Matrix<T>&);
     void operator+= (const Matrix<T>&);
     void operator-= (const Matrix<T>&);
     void operator*= (T);
     void operator+= (T);
     void operator-= (T);
     */
    
};

//==========================================================CONSTRUCTORS============================================================
template <typename T>
inline Matrix<T>::Matrix()
{
    this->_rows = _DEFAULT_MATRIX_SIZE;
    this->_cols = _DEFAULT_MATRIX_SIZE;
    _initilalizeMatrix();
}

template<typename T>
inline Matrix<T>::Matrix(const std::string& file_path)
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
    
    file_reader >> this->_rows;
    file_reader >> this->_cols;
    
    _initilalizeMatrix();
    file_reader >> tmp_x;      //skip the line number
    
    while (!file_reader.eof())
    {
        file_reader >> tmp_x;
        file_reader >> tmp_y;
        (*this)(tmp_x - 1 ,tmp_y - 1, 1);
    }
    
    file_reader.close();
}

template <typename T>
inline Matrix<T>::Matrix(int rows, int cols)
{
    this->_rows = rows;
    this->_cols = cols;
    _initilalizeMatrix();
}

template <typename T>
inline Matrix<T>::Matrix(const Matrix<T>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
template <typename T>
inline Matrix<T>::~Matrix()
{
    if (this->_edges != NULL)
    {
        delete [] this->_edges;
    }
    
    if (sparse_form != NULL)
    {
        for(int i=0; i < sparse_form_size; i++)
        {
            delete sparse_form[i];
        }
        delete [] sparse_form;
    }
}

//===========================================================ACCESSORS===============================================================
template <typename T>
inline int Matrix<T>::getNumberOfRows()
{
    return this->_rows;
}

template <typename T>
inline int Matrix<T>::getNumberOfColumns()
{
    return this->_cols;
}

template <typename T>
inline bool Matrix<T>::isSquare() const
{
    return (this->_rows == this->_cols);
}

template <typename T>
inline bool Matrix<T>::isSymmetric() const
{
    //sym. matrix has to be square
    if (this->_rows != this->_cols)
    {
        return false;
    }
    
    //chaking for entries to be equal
    for(int i = 0; i < this->_rows; i++)
    {
        for(int j = 0; j < this->_cols; j++)
        {
            if (this->_edges[(i * this->_cols) + j] != this->_edges[(j * this->_cols) + i])
            {
                return false;
            }
        }
    }
    
    return true;
}

template <typename T>
inline sparse_matrix_element<T>** Matrix<T>::getSparseForm()
{
    if(sparse_form != NULL)
    {
        return this->sparse_form;
    }
    
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        if(this->_edges[i] != 0)
        {
            sparse_form_size++;
        }
    }
    
    sparse_form = new sparse_matrix_element<T>* [sparse_form_size];
    int counter=0;
    
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        if(this->_edges[i] != 0)
        {
            sparse_matrix_element<T> *ele= new sparse_matrix_element<T>;
            ele->row_index = i / this->_rows;
            ele->col_index = i % this->_rows;
            ele->value = this->_edges[i];
            sparse_form[counter]=ele;
            counter++;
        }
    }
    
    return this->sparse_form;
}

template <typename T>
inline Matrix<T> Matrix<T>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B)
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
    //Initializing and allocating the product matrix
    Matrix<T> res_matrix(num_in_A, num_in_B);
    
    int counter = 0;
    
    for (int i=0; i < vec_A.size(); i++)
    {
        for(int j=0; j < vec_B.size(); j++)
        {
            if ( vec_A[i] == 1 && vec_B[j] ==1)
            {
                res_matrix._edges[counter] = (*this)(i,j);
                counter++;
            }
        }
    }
    return res_matrix;
}

template <typename T>
inline std::vector<int> Matrix<T>::getNeighbors(int vertex)
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
inline std::vector<T> Matrix<T>::getTopEigenVector()
{
#ifdef __linux__
    int arr_size = (0.5 * this->_rows * (this->_rows+1));
    T sym_edges[arr_size];
    int counter = 0;
    
    for(int i = 0; i < this->_rows; i++)
    {
        for(int j = i; j < this->_cols; j++)
        {
            sym_edges[counter] = (*this)(i,j);
            counter++;
        }
    }
    
    ARdsSymMatrix<T> ARMatrix(this->_rows, sym_edges, 'L');
    ARluSymStdEig<T> eigProb(1, ARMatrix, "LM", 10);
    eigProb.FindEigenvectors();
    
    std::vector<T> eigen_vec (eigProb.GetN());
    for (int i=0; i < eigProb.GetN() ; i++)
    {
        eigen_vec[i] = eigProb.Eigenvector(0 ,i);
    }
    
    return eigen_vec;
#endif
}

template <typename T>
inline int Matrix<T>:: getSparseFormSize()
{
    if(sparse_form == NULL)
    {
        this->getSparseForm();
    }
    return this->sparse_form_size;
}

/*
 * returns a vector of type T where each entry corresponds to sum of entries in a matrix row.
 * Delete this object after usage.
 * @return std::vector<T>
 */
template <typename T>
inline std::vector<T> Matrix<T>::getSumOfRows()
{
    std::vector<T> sum_vector(this->_rows);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        sum_vector[(i / this->_rows)] += this->_edges[i];
    }
    
    return sum_vector;
}

//===========================================================MUTATORS================================================================



template <typename T>
inline void Matrix<T>::insert(int i, int j, T value)
{
    if ( (i >= _rows) || (i < 0) || (j >= _cols) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    
    this->_edges[(i * this->_cols) + j] = value;
}

//==========================================================OPERATORS================================================================
template <typename T>
inline Matrix<T> Matrix<T>::kron(const Matrix<T>& matrix)
{
    // checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    //Initializing and allocating the product matrix
    int prod_size = this->_rows * matrix._rows;
    Matrix<T> prod_matrix(prod_size, prod_size);
    
    /*
     *  Calculating the kronecker product:
     *  The indices of the product matrix is calculated by:
     *      i = (i_outer*size) + i_inner
     *      j = (j_outer*size) + j_inner
     */
    for (int i_outer = 0; i_outer < this->_rows; i_outer++)
    {
        for (int j_outer=0; j_outer < this->_cols; j_outer++)
        {
            for(int i_inner=0; i_inner < matrix._rows; i_inner++)
            {
                for(int j_inner=0; j_inner < matrix._cols; j_inner++)
                {
                    prod_matrix.insert((i_outer*matrix._rows) + i_inner,
                                       (j_outer*matrix._cols) + j_inner,
                                       matrix(i_inner, j_inner) * (*this)(i_outer , j_outer));
                }
            }
        }
    }
    
    //new method to take more advantage of 1d array and possibly faster, the assignment needs some thinking
    
    //    Matrix<T>* kron_prod = new Matrix<T>(prod_size, prod_size);
    //    int counter = 0;
    //    for (int i = 0; i < this->_getArrSize(); i++)
    //    {
    //        for (int j = 0; j < matrix._getArrSize(); j++)
    //        {
    //            kron_prod->_edges[counter++] = this->_edges[i] * matrix._edges[j];
    //        }
    //    }
    //    std::cout << (*kron_prod) << std::endl;
    //    std::cout << ((*kron_prod) == (*prod_matrix)) << std::endl;
    
    return prod_matrix;
}

template <typename T>
inline Matrix<T> Matrix<T>::diagonalVectorTimesMatrix(const std::vector<T>& vec)
{
    if(this->_rows != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    Matrix<T> ret_matrix(*this);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        ret_matrix._edges[i] = vec[i/this->_rows] * this->_edges[i];
    }
    
    return ret_matrix;
}

template <typename T>
inline Matrix<T> Matrix<T>::matrixTimesDiagonalVector(const std::vector<T>& vec)
{
    if(this->_cols != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    Matrix<T> ret_matrix(*this);
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

template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix)
{
    stream<< "Size: " << matrix._rows << "*" << matrix._cols << '\n';
    for (int i = 1; i <= matrix._getArrSize(); i++)
    {
        stream << matrix._edges[i-1] << ' ';
        if ( (i % matrix._rows) == 0)
        {
            stream << "\n";
        }
    }
    stream << "\n\n\n";
    return stream;
}
    
template <typename T>
inline T Matrix<T>::operator()(int i, int j) const
{
    if ( (i >= _rows) || (i < 0) || (j >= _cols) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    
    return this->_edges[(i * this->_cols) + j];
}
    
template <typename T>
inline void Matrix<T>::operator()(int i, int j, T value)
{
    if ( (i >= _rows) || (i < 0) || (j >= _cols) || (j < 0) )
    {
        throw IndexOutOfBoundsException();
    }
    
    this->_edges[(i * this->_cols) + j] = value;
}
    
template <typename T>
inline Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrix)
{
    if (this->_edges != NULL)
    {
        delete [] this->_edges;
    }
    _copy(matrix);
}

template <typename T>
inline bool Matrix<T>::operator==(const Matrix<T>& matrix)
{
    // checking for dimension equality
    if ((this->_cols != matrix._cols) || (this->_rows != matrix._rows))
    {
        return false;
    }
    
    //checking for entry equality
    for (int i = 0; i < this->_getArrSize(); i++)
    {
        if ( this->_edges[i] != matrix._edges[i])
        {
            return false;
        }
    }
    
    return true;
}
    
//===========================================================PRIVATE=================================================================
template <typename T>
inline void Matrix<T>::_copy(const Matrix<T>& matrix)
{
    this->_rows = matrix._rows;
    this->_cols = matrix._cols;
    _initilalizeMatrix();
    
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        this->_edges[i] = matrix._edges[i];
    }
}
    
template <typename T>
inline void Matrix<T>::_initilalizeMatrix()
{
    this->sparse_form_size = 0;
    this->sparse_form = NULL;
    this->_edges = new T[this->_getArrSize()];
    
    //checking for memory issues.
    if (this->_edges == NULL)
    {
        throw OutOfMemoryException();
    }
    
    for ( int i = 0; i < this->_getArrSize(); i++)
    {
        this->_edges[i] = 0;
    }
}
    
template<typename T>
inline size_t Matrix<T>::_getArrSize() const
{
    return (this->_rows * this->_cols);
}
    
//===================================================================================================================================
    
    
#endif
    