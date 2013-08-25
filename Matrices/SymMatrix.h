/************************************************************************************
 * Symmetrix Matrix Data Structure. This structure uses an array to store the values* 
 * of a lower triagular matrix.                                                     *
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

template <typename DT>
class SymMatrix;

template <typename DT>
std::ostream& operator<< (std::ostream&, const SymMatrix<DT>&);

/*
 * SymMatrix class definition and method declarations.
 */
template <typename DT>
class SymMatrix
{
private:
    /*
     * Private methods (used internally)
     */
    size_t _getArrSize() const;
    void _initializeMatrix(bool);
    void _copy(const SymMatrix<DT>&);

    /*
     * Private constants.
     */
    static const int _DEFAULT_MATRIX_SIZE;
    static const  T _DEFAULT_MATRIX_ENTRY;
    
protected:
    int _size;
    DT* _edges;

public:
      
    /**************
     *Constructors*
     **************/
    SymMatrix(bool fill = true);
    SymMatrix(int size, bool fill = true);
    SymMatrix(const SymMatrix<DT>&);
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
    int getNumberOfRows() const;
    int getNumberOfColumns() const;
    std::vector<SparseElement<T> >getSparseForm() const;
    DenseMatrix1D<DT> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B);
    
    /**********
     *MUTATORS*
     **********/
    
    /**********
    *OPERATIONS*
    **********/
    T getFrobNorm() const;
    T* getTopEigenVector() const;
    SymMatrix<T> transpose() const;
    std::vector<T> getSumOfRows() const;
    SymMatrix<T> kron(const SymMatrix<T>& matrix) const;
    SymMatrix<T> diagonalVectorTimesMatrix(const std::vector<T>&) const;
    SymMatrix<T> matrixTimesDiagonalVector(const std::vector<T>&) const;

    /*************
    *  MPI Send  *
    **************/
    #ifdef USE_MPI
    void MPI_Send_Matrix(int, int);
    void MPI_Bcast_Send_Matrix(int);
    #endif

    /**********
     *OPERATORS*
     **********/
    T& operator()(int, int);
    void operator=(const SymMatrix<T>&);
    bool operator==(const SymMatrix<T>&); 
    SymMatrix<T> operator+(const SymMatrix<T>& other_matrix) const;
    SymMatrix<T> operator-(const SymMatrix<T>& other_matrix) const;
    SymMatrix<T> operator*(const SymMatrix<T>& other_matrix) const;
    friend std::ostream& operator<< <> (std::ostream& stream, const SymMatrix<T>& matrix);
};

//==========================================================CONSTANTS============================================================
template <typename T>
const int SymMatrix<T>::_DEFAULT_MATRIX_SIZE = 1;
template <typename T>
const T SymMatrix<T>::_DEFAULT_MATRIX_ENTRY = 1;
//==========================================================CONSTRUCTORS============================================================

/*
 * Default constructor:
 * Construct a matrix of size _DEFAULT_MATRIX_SIZE * _DEFAULT_MATRIX_SIZE initialized to 0.
 * @pram bool fill: fills the matrix with 0's. default value is true
 */
template <typename DT>
inline SymMatrix<DT>::SymMatrix(bool fill = true)
{
    this->_size = _DEFAULT_MATRIX_SIZE;
    _initializeMatrix(fill);
}

/*
 * SymMatrix constructor:
 * Construct a matrix by reading a matrix file, specified in readme.txt file the nodes that exist will have _DEFAULT_MATRIX_ENTRY value.
 * @pram std::string : path to the file
 */
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
        throw FileDoesNotExistException(file_path);
    }
    
    file_reader >> this->_size;
    file_reader >> tmp_x;
    
    
    if (this->_size != tmp_x)
    {
        file_reader.close();
        throw NotASquareMatrixException();
    }
    
    _initializeMatrix(false);
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
template <typename DT>
inline SymMatrix<DT>::SymMatrix(int size, bool fill = true)
{
    this->_size = size;
    _initializeMatrix(fill);
}

#ifdef USE_MPI
/*
 * SymMatrix constructor:
 * Construct a matrix by receiving matrix data from a MPI_Send_Matrix call.
 * Note the matrix being send should be of type DenseMatrix1D, DenseMatrix2D, SymmMatrix
 * @pram int source: Sender's ID
 * @pram int tag: sender's tag
 * @pram MPI_Status: MPI_Status object
 */
template <typename DT>
inline SymMatrix<DT>::SymMatrix(int source, int tag, MPI_Status& stat)
{
    MPI_Recv(&this->_size, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
    _initializeMatrix(false);
    MPI_Recv(this->_edges, this->_getArrSize()*sizeof(DT), MPI_BYTE, source, tag + 2, MPI_COMM_WORLD, &stat);
}

/*
 * SymMatrix constructor:
 * Construct a matrix by receiving matrix data from a MPI_Bcast_Send_Matrix call.
 * Note the matrix being send should be of type DenseMatrix1D, DenseMatrix2D, SymmMatrix
 * @pram int source: Sender's ID
 * @pram MPI_Status: MPI_Status object
 */
template <typename DT>
inline SymMatrix<DT>::SymMatrix(int source, MPI_Status& stat)
{
    MPI_Bcast(&this->_size, 1, MPI_INT, source, MPI_COMM_WORLD);
    _initializeMatrix(false);
    MPI_Bcast(this->_edges, this->_getArrSize()*sizeof(DT), MPI_BYTE, source, MPI_COMM_WORLD);
}

#endif

/*
 * SymMatrix copy constructor:
 * Construct a copy of a matrix.
 * @pram DenseMatrix2D<T>
 */
template <typename DT>
inline SymMatrix<DT>::SymMatrix(const SymMatrix<DT>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
/*
 * SymMatrix destructor:
 * deleted the arrays used in the matrix. NOTE: user is responsible for deleting the objects in the matrix if they are dynamically allocated.
 */
template <typename DT>
inline SymMatrix<DT>::~SymMatrix()
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
inline std::vector<SparseElement<T> > DenseMatrix1D<T>::getSparseForm() const
{    
    std::vector<SparseElement<T> > sparse_form;
    for(int i = 0; i < this->_rows; i++)
    {
        for(int j = 0; j < this->_cols; j++)
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
template <typename DT>
inline DenseMatrix1D<DT> SymMatrix<DT>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B)
{
    // int num_in_A = 0;
    // for (int i=0; i< vec_A.size(); i++)
    // {
    //     if (vec_A[i] != 0)
    //     {
    //         num_in_A++;
    //     }
    // }
    // int num_in_B = 0;
    // for (int i=0; i < vec_B.size(); i++)
    // {
    //     if( vec_B[i] != 0)
    //     {
    //         num_in_B++;
    //     }
    // }
    // //Initializing and allocating the product matrix
    // Matrix<DT> res_matrix(num_in_A, num_in_B);
    // int counter = 0;
    // for (int i=0; i< vec_A.size(); i++)
    // {
    //     for(int j=0; j< vec_B.size(); j++)
    //     {
    //         if ( vec_A[i] == 1 && vec_B[j] ==1)
    //         {
    //             res_matrix.insert(counter/num_in_B, counter%num_in_B, (*this)(i, j));
    //             counter++;
    //         }
    //     }
    // }
    // return res_matrix;
}
//===========================================================MUTATORS==================================================================

//===========================================================OPERATIONS================================================================
template <typename DT> //should be a better implementation
inline SymMatrix<DT> SymMatrix<DT>::kron(const SymMatrix<DT>& matrix)
{

    
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
        for (int j_outer = i_outer; j_outer < this->_size; j_outer++)
        {
            for(int i_inner = 0; i_inner < matrix._size; i_inner++)
            {
                for(int j_inner = 0; j_inner < matrix._size; j_inner++)
                {
                    prod_matrix((i_outer*matrix._size) + i_inner,
                                        (j_outer*matrix._size)+j_inner) = 
                                        (*this)(i_outer, j_outer) * matrix(i_inner, j_inner);
                }
            }
        }
    }
    
    return prod_matrix;
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
#ifdef ARPACK
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

//===========================================================MPI SEND/REC================================================================
#ifdef USE_MPI
template <typename DT>
void SymMatrix<DT>::MPI_Send_SymMatrix(int dest, int tag)
{
    MPI_Send(&this->_size, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
    MPI_Send(this->_edges, this->_getArrSize()*sizeof(DT), MPI_BYTE, dest, tag + 2, MPI_COMM_WORLD);
}

template <typename DT>
void SymMatrix<DT>::MPI_Bcast_Send_SymMatrix(int source)
{
    MPI_Bcast(&this->_size, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(this->_edges, this->_getArrSize()*sizeof(DT), MPI_BYTE, source, MPI_COMM_WORLD);
}
#endif
//==========================================================OPERATIONS================================================================
template <typename DT>
inline SymMatrix<DT> SymMatrix<DT>::diagonalVectorTimesMatrix(const std::vector<DT>& vec)
{
    if(this->_rows != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    SymMatrix<DT> ret_matrix(*this);
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
    
    SymMatrix<DT> ret_matrix(*this);
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
inline std::ostream& operator<<(std::ostream& stream, const SymMatrix<DT>& matrix)
{
    stream<< "Size: " << matrix.getSize() << "*" << matrix.getSize()<< '\n';
    for (int i = 0; i < matrix.getSize(); i++)
    {
        for( int j = 0; j< matrix.getSize(); j++)
        {
            stream << matrix(i, j) << ' ';
        }
        stream << "\n";
    }
    stream << "\n\n\n";
    return stream;
}

template <typename DT>
inline DT& SymMatrix<DT>::operator()(int i, int j) const
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
    _initializeMatrix();
  
    for(int i=0; i < this->_getArrSize() ; i++)
    {
        this->_edges[i] = matrix._edges[i];
    }
}

template <typename DT>
inline void SymMatrix<DT>::_initializeMatrix()
{
    try
    {
        this->_edges = new DT[this->_getArrSize()];    
    } 
    catch(std::bad_alloc& e) 
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
