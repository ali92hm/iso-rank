/************************************************************************************
 * SymMatrix Matrix Data Structure. This structure uses an array to store the values* 
 * of a lower triangular matrix.                                                     *
 ************************************************************************************/

#ifndef _SymMatrix_h
#define _SymMatrix_h


#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "MatrixExceptions.h"
#include "SparseElement.h"
#include "DenseMatrix1D.h"

#ifdef EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#ifdef ARPACK
#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

template <typename T>
class SymMatrix;

template <typename T>
std::ostream& operator<< (std::ostream&, const SymMatrix<T>&);

/*
 * SymMatrix class definition and method declarations.
 */
template <typename T>
class SymMatrix
{
private:
    /*
     * Private methods (used internally)
     */
    size_t _getArrSize() const;
    void _initializeMatrix(bool);
    void _copy(const SymMatrix<T>&);

    /*
     * Private constants.
     */
    static const int _DEFAULT_MATRIX_SIZE;
    static const  T _DEFAULT_MATRIX_ENTRY;
#ifdef USE_MPI
    static const int _SPARSE_FORM;
    static const int _DENSE_FORM;
    static const int _SYM_SPARSE_FORM;
    static const int _SYM_DENSE_FORM;
#endif
    
protected:
    int _size;
    T* _edges;

public:
      
    /**************
     *Constructors*
     **************/
    SymMatrix(bool fill = true);
    SymMatrix(int size, bool fill = true);
    SymMatrix(const SymMatrix<T>&);
    SymMatrix(const std::string&);

    #ifdef USE_MPI
    SymMatrix(int,MPI_Status&);
    SymMatrix(int,int,MPI_Status&);
    #endif

    /************
     *Destructor*
     ************/
    virtual ~SymMatrix();
    
    /***********
     *ACCESSORS*
     ***********/
    int getNumberOfRows();
    int getNumberOfColumns();
    std::vector<SparseElement<T> >getSparseForm();
    DenseMatrix1D<T> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int>& vec_B);
    
    /**********
     *MUTATORS*
     **********/
    
    /**********
    *OPERATIONS*
    **********/
    T getFrobNorm();
    T* getTopEigenVector();
    SymMatrix<T> transpose();
    std::vector<T> getSumOfRows();
    std::vector<int> getNeighbors(int vertex);
    SymMatrix<T> kron(const SymMatrix<T>& matrix);
    DenseMatrix1D<T> diagonalVectorTimesMatrix(const std::vector<T>&);
    DenseMatrix1D<T> matrixTimesDiagonalVector(const std::vector<T>&);

    /*************
    *  MPI Send  *
    **************/
    #ifdef USE_MPI
    void MPI_Send_Matrix(int, int, bool sparse = false);
    void MPI_Bcast_Send_Matrix(int, bool sparse = false);
    #endif

    /**********
     *OPERATORS*
     **********/
    T& operator()(int, int);
    void operator=(const SymMatrix<T>&);
    bool operator==(const SymMatrix<T>&); 
    DenseMatrix1D<T> operator+(const DenseMatrix1D<T>& other_matrix);
    DenseMatrix1D<T> operator-(const DenseMatrix1D<T>& other_matrix);
    DenseMatrix1D<T> operator*(const DenseMatrix1D<T>& other_matrix);
    friend std::ostream& operator<< <> (std::ostream& stream, const SymMatrix<T>& matrix);
};

//==========================================================CONSTANTS============================================================
template <typename T>
const int SymMatrix<T>::_DEFAULT_MATRIX_SIZE = 1;
template <typename T>
const T SymMatrix<T>::_DEFAULT_MATRIX_ENTRY = 1;
#ifdef USE_MPI
template <typename T>
const int SymMatrix<T>::_SENDING_DENSE_FORM = 0;
template <typename T>
const int SymMatrix<T>::_SENDING_SPARSE_FORM = 1;
template <typename T>
const int SymMatrix<T>::_SENDING_SYM_DENSE_FORM = 2;
template <typename T>
const int SymMatrix<T>::_SENDING_SYM_SPARSE_FORM = 3;
#endif
//==========================================================CONSTRUCTORS============================================================

/*
 * Default constructor:
 * Construct a matrix of size _DEFAULT_MATRIX_SIZE * _DEFAULT_MATRIX_SIZE initialized to 0.
 * @pram bool fill: fills the matrix with 0's. default value is true
 */
template <typename T>
inline SymMatrix<T>::SymMatrix(bool fill)
{
    this->_size = _DEFAULT_MATRIX_SIZE;
    _initializeMatrix(fill);
}

/*
 * SymMatrix constructor:
 * Construct a matrix by reading a matrix file, specified in readme.txt file the nodes that exist will have _DEFAULT_MATRIX_ENTRY value.
 * @pram std::string : path to the file
 */
template<typename T>
inline SymMatrix<T>::SymMatrix(const std::string& file_path)
{
    int tmp_x;
    int tmp_y;
    std::ifstream file_reader;
    file_reader.open(file_path.c_str());
    
    if(file_reader.fail())
    {
        file_reader.close();
        throw FileDoesNotExistException(file_path);
    }
    
    file_reader >> this->_size;
    file_reader >> tmp_x;
    
    
    if (this->_size != tmp_x)
    {
        file_reader.close();
        throw NotASquareMatrixException();
    }
    
    _initializeMatrix(true);
    file_reader >> tmp_x;      //skip the line number
    
    while (!file_reader.eof())
    {
        file_reader >> tmp_x;
        file_reader >> tmp_y;
        (*this)(tmp_x - 1 ,tmp_y - 1) = 1;
    }
    
    file_reader.close();
}

/*
 * SymMatrix constructor:
 * Construct a matrix n*n matrix initialized to 0.
 * @pram int size: number of rows and columns
 * @pram bool fill: fills the matrix with 0's. default value is true
 */
template <typename T>
inline SymMatrix<T>::SymMatrix(int size, bool fill)
{
    this->_size = size;
    _initializeMatrix(fill);
}

#ifdef USE_MPI
/*
 * SymMatrix constructor:
 * Construct a matrix by receiving matrix data from a MPI_Send_Matrix call.
 * Note the matrix being send should be of type DenseMatrix1D, SymMatrix, SymmMatrix
 * @pram int source: Sender's ID
 * @pram int tag: sender's tag
 * @pram MPI_Status: MPI_Status object
 */
template <typename T>
inline SymMatrix<T>::SymMatrix(int source, int tag, MPI_Status& stat)
{
    MPI_Recv(&this->_size, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
    _initializeMatrix(false);
    MPI_Recv(this->_edges, this->_getArrSize()*sizeof(T), MPI_BYTE, source, tag + 2, MPI_COMM_WORLD, &stat);
}

/*
 * SymMatrix constructor:
 * Construct a matrix by receiving matrix data from a MPI_Bcast_Send_Matrix call.
 * Note the matrix being send should be of type DenseMatrix1D, SymMatrix, SymmMatrix
 * @pram int source: Sender's ID
 * @pram MPI_Status: MPI_Status object
 */
template <typename T>
inline SymMatrix<T>::SymMatrix(int source, MPI_Status& stat)
{
    MPI_Bcast(&this->_size, 1, MPI_INT, source, MPI_COMM_WORLD);
    _initializeMatrix(false);
    MPI_Bcast(this->_edges, this->_getArrSize()*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD);
}

#endif

/*
 * SymMatrix copy constructor:
 * Construct a copy of a matrix.
 * @pram SymMatrix<T>
 */
template <typename T>
inline SymMatrix<T>::SymMatrix(const SymMatrix<T>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
/*
 * SymMatrix destructor:
 * deleted the arrays used in the matrix. NOTE: user is responsible for deleting the objects in the matrix if they are dynamically allocated.
 */
template <typename T>
inline SymMatrix<T>::~SymMatrix()
{
    if (_edges != NULL)
    {
        delete [] _edges;
    }
}

//===========================================================ACCESSORS===============================================================
/*
 * Returns the number of the rows.
 */
template <typename T>
inline int SymMatrix<T>::getNumberOfRows()
{
    return this->_size;
}

/*
 * Returns the number of the columns.
 */
template <typename T>
inline int SymMatrix<T>::getNumberOfColumns()
{
    return this->_size;
}

/*
 * Returns a std::vector<SparseElement<T>> of SparseElement objects that contain the i, j and value of the none-zero edges.
 */
template <typename T>
inline std::vector<SparseElement<T> > SymMatrix<T>::getSparseForm()
{    
    std::vector<SparseElement<T> > sparse_form;
    for(int i = 0; i < this->_size; i++)
    {
        for(int j = 0; j < this->_size; j++)
        {
            if((*this)(i, j) != 0)
            {
                sparse_form.push_back(SparseElement<T>(i,j, (*this)(i, j)));
            }
        }
    }
    return sparse_form;
}

/*
 * Returns a DenseMatrix1D that has the rows that are marked 1 in vec_A, columns that are marked 1 in vec_B
 * @pram std::vector<int>: vector of 0's and 1's for row selection
 * @pram std::vector<int>: vector of 0's and 1's for column selection
 */
template <typename T>
inline DenseMatrix1D<T> SymMatrix<T>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int>& vec_B)
{
}
//===========================================================MUTATORS==================================================================

//===========================================================OPERATIONS================================================================
/*
 * Returns the frobenius norm.
 */
template <typename T>
inline T SymMatrix<T>::getFrobNorm()
{

}

/*
 * Returns a std::vector<T> that contains the values of the eigenvector associated to the largest eigenvalue
 */
template <typename T>
inline T* SymMatrix<T>::getTopEigenVector()
{
#ifdef ARPACK
    ARdsSymMatrix<T> ARMatrix(this->_size, this->_edges, 'L');
    ARluSymStdEig<T> eigProb(1, ARMatrix, "LM", 10);
    eigProb.FindEigenvectors();
    T* eigen_vector = new T[this->_size];
     
    for (int i=0; i < this->_size; i++)
    {
        eigen_vector[i] = eigProb.Eigenvector(0,i);
    }
    
    return eigen_vector;
#endif

#ifdef EIGEN

#endif
}

/*
 * Returns a SymMatrix<T> of the transpose of this (creates a copy)
 */
template <typename T>
inline SymMatrix<T> SymMatrix<T>::transpose()
{
}

/*
 * Returns a std::vector<T> that contains sum of the values in each row
 */
template <typename T>
inline std::vector<T> SymMatrix<T>::getSumOfRows()
{
}

/*
 * Returns a std::vector<int> of the elements with value of 1 in a columns.
 * @pram int vertex
 */
template <typename T>
inline std::vector<int> SymMatrix<T>::getNeighbors(int vertex)
{
}

/*
 * Returns a SymMatrix that is the kronecker product of this and another matrix.
 * @pram SymMatrix 
 */
template <typename T> 
inline SymMatrix<T> SymMatrix<T>::kron(const SymMatrix<T>& matrix)
{
}

/*
 * Returns a SymMatrix<T> of a vector that contains the diagonal values of a diagonal matrix times this matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename T>
inline DenseMatrix1D<T> SymMatrix<T>::diagonalVectorTimesMatrix(const std::vector<T>& vec)
{
}

/*
 * Returns a SymMatrix<T> of the product of this matrix a vector that contains the diagonal values of a diagonal matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename T>
inline DenseMatrix1D<T> SymMatrix<T>::matrixTimesDiagonalVector(const std::vector<T>& vec)
{
}


//===========================================================MPI SEND/REC================================================================
#ifdef USE_MPI
/*
 * Sends a SymMatrix<T> using MPI_Send.
 * @pram int destination processor
 * @pram int tag 
 * @pram bool sending sparse form, default value is false.
 */
template <typename T>
void SymMatrix<T>::MPI_Send_SymMatrix(int dest, int tag, bool sparse = false)
{
    if(sparse)
    {
        MPI_Send(&_SYM_SPARSE_FORM, 1 ,MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        std::vector<SparseElement<T> > sparse_form;
        for (int i = 0; i < this->_size; i++)
        {
            for (int j = i ; j < this->_size; j++)
            {   
                sparse_form.push_back(SparseElement(i,j,(*this)(i,j));
            }
        }
        MPI_Send(&sparse_form.size(), 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<T>), MPI_BYTE, dest, tag + 3, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&_SYM_DENSE_FORM, 1 ,MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        MPI_Send(&this->_size, 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(this->_edges, this->_getArrSize()*sizeof(T), MPI_BYTE, dest, tag + 3, MPI_COMM_WORLD);
    }
}

/*
 * Sends a SymMatrix<T> using MPI_Bcast.
 * @pram int source sender's ID
 * @pram bool sending sparse form, default value is false.
 */
template <typename T>
void SymMatrix<T>::MPI_Bcast_Send_SymMatrix(int source, bool sparse = false)
{
    if(sparse)
    {
        MPI_Bcast(&_SYM_SPARSE_FORM, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        std::vector<SparseElement<T> > sparse_form;
        for (int i = 0; i < this->_size; i++)
        {
            for (int j = i ; j < this->_size; j++)
            {   
                sparse_form.push_back(SparseElement(i,j,(*this)(i,j));
            }
        }
        MPI_Bcast(&sparse_form.size(), 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<T>), MPI_BYTE, source, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&_SYM_DENSE_FORM, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_size, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(this->_edges, this->_getArrSize()*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD);
    }

}
#endif
//==========================================================OPERATIONS================================================================


//==========================================================OPERATORS================================================================
/*
 * Overloaded () operator for accessing and changing the values inside a matrix
 * @pram: int i
 * @pram: int j
 */
template <typename T>
inline T& SymMatrix<T>::operator()(int i, int j)
{
    if(i<=j)
    {
        return this->_edges[i+size_t(j)*(j+1)/2];
    }
    else
    {
        return this->_edges[j+size_t(i)*(i+1)/2];
    }
}

/*
 * Overloaded = operator copies the content of another matrix to this
 * @pram: SymMatrix<T> 
 */
template <typename T>
inline void SymMatrix<T>::operator=(const SymMatrix<T>& matrix)
{
    delete [] this->_edges;
    _copy(matrix);
}

/*
 * Overloaded == operator to compare two matrices
 * @pram: SymMatrix<T> 
 */
 template <typename T>
 inline bool SymMatrix<T>::operator==(const SymMatrix<T>& matrix)
 {
    if (this->_size != matrix._size)
    {
        return false;
    }

    for(int i = 0; i < this->_getArrSize(); i++)
    {
        if (this->_edges[i] != matrix._edges[i])
        {
            return false;
        }
    }
    return true;
 }

/*
 * Returns a DenseMatrix1D that is the sum of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> SymMatrix<T>::operator+(const DenseMatrix1D<T>& other_matrix)
{
}

/*
 * Returns a DenseMatrix1D that is the difference of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> SymMatrix<T>::operator-(const DenseMatrix1D<T>& other_matrix)
{
}

/*
 * Returns a DenseMatrix1D that is product of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> SymMatrix<T>::operator*(const DenseMatrix1D<T>& other_matrix)
{
}

/*
 * Overloaded ostream operator for printing a matrix
 * @pram: std::ostream 
 * @pram: SymMatrix<T>
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const SymMatrix<T>& matrix)
{
    stream<< "Size: " << matrix._size << "*" << matrix._size<< '\n';
    for (int i = 0; i < matrix._size; i++)
    {
        for( int j = 0; j< matrix._size; j++)
        {
            if(i<=j)
            {
                stream <<  matrix._edges[i+size_t(j)*(j+1)/2] << ' ';;
            }
            else
            {
                stream <<  matrix._edges[j+size_t(i)*(i+1)/2] << ' ';;
            }
        }
        stream << "\n";
    }
    stream << "\n\n\n";
    return stream;
}

//===========================================================PRIVATE=================================================================
/*
 * Make a deep copy of a matrix object
 * @pram: SymMatrix<T> 
 */
template <typename T>
inline void SymMatrix<T>::_copy(const SymMatrix<T>& matrix)
{
    this->_size = matrix._size;
    _initializeMatrix(false);
    memcpy(this->_edges, matrix._edges, this->_getArrSize() * sizeof(T));
}

/*
 * makes a 1D array of size _rows*_cols
 * @pram: bool fill: if true: initialize values to 0.
 */  
template <typename T>
inline void SymMatrix<T>::_initializeMatrix(bool fill)
{
    try
    {
        if (fill)
        {
            this->_edges = new T[this->_getArrSize()]();
        }
        else
        {
            this->_edges = new T[this->_getArrSize()];
        }
    }
    catch (std::bad_alloc& e)
    {
        throw OutOfMemoryException();
    }
}

/*
 * Return the size of the internal array
 * @return: size_t
 */
template<typename T>
inline size_t SymMatrix<T>::_getArrSize() const
{
    return (int)(0.5 * this->_size * (this->_size + 1));
}

//===================================================================================================================================

#endif
