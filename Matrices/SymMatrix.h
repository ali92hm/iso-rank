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
#ifdef USE_MPI
    static const int SENDING_SPARSE_FORM;
    static const int SENDING_DENSE_FORM:
    static const int SENDING_SYM_SPARSE_FORM;
    static const int SENDING_SYM_DENSE_FORM:
#endif
    
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
    void MPI_Send_Matrix(int, int, bool);
    void MPI_Bcast_Send_Matrix(int, bool);
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
    DenseMatrix1D<DT> res_matrix(num_in_A, num_in_B);
    int counter = 0;
    for (int i=0; i< vec_A.size(); i++)
    {
        for(int j=0; j< vec_B.size(); j++)
        {
            if ( vec_A[i] == 1 && vec_B[j] ==1)
            {
                res_matrix(counter/num_in_B, counter%num_in_B) = (*this)(i, j);
                counter++;
            }
        }
    }
    return res_matrix;
}
//===========================================================MUTATORS==================================================================

//===========================================================OPERATIONS================================================================
/*
 * Returns the frobenius norm.
 */
template <typename T>
inline T SymMatrix<T>::getFrobNorm() const
{

}

/*
 * Returns a std::vector<T> that contains the values of the eigenvector associated to the largest eigenvalue
 */
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

#ifdef EIGEN
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> A_eigen = Eigen::MatrixXd::Zero(this->_size, this->_size);
        for (int i=0; i< this->_size; i++) 
        {
            for(int j = 0; j < this->_size; j++)
            {
                A_eigen(i,j) = (*this)(i,j);
            }
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        es.compute(A_eigen);
        
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> evals_eigen = es.eigenvalues();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> evecs_eigen = es.eigenvectors();
        T* eigen_vector = new T[this->_size];

        for ( int i=0; i < this->_size; i++)
        {   
            if (evals_eigen(i) == evals_eigen.maxCoeff())
            {
                for (int j= 0; j < this->_size ; j++)
                {
                    eigen_vector[j] = evecs_eigen(j,i);
                }
                break;
            }
        }
    return eigen_vector;
#endif
}

/*
 * Returns a SymMatrix<T> of the transpose of this (creates a copy)
 */
template <typename T>
inline SymMatrix<T>::transpose() const
{
    return SymMatrix<T>(*this);
}

/*
 * Returns a std::vector<T> that contains sum of the values in each row
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

/*
 * Returns a std::vector<int> of the elements with value of 1 in a columns.
 * @pram int vertex
 */
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

/*
 * Returns a DenseMatrix2D that is the kronecker product of this and another matrix.
 * @pram DenseMatrix2D 
 */
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

/*
 * Returns a DenseMatrix2D<T> of a vector that contains the diagonal values of a diagonal matrix times this matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename DT>
inline SymMatrix<DT> SymMatrix<DT>::diagonalVectorTimesMatrix(const std::vector<DT>& vec)
{
    // if(this->_rows != vec.size())
    // {
    //     throw DimensionMismatchException();
    // }
    
    // SymMatrix<DT> ret_matrix(*this);
    // for(int i = 0; i < this->_getArrSize(); i++)
    // {
    //     ret_matrix._edges[i] = vec[i/this->_rows] * this->_edges[i];
    // }
    
    // return ret_matrix;
}

/*
 * Returns a DenseMatrix2D<T> of the product of this matrix a vector that contains the diagonal values of a diagonal matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename DT>
inline SymMatrix<DT> SymMatrix<DT>::matrixTimesDiagonalVector(const std::vector<DT>& vec)
{
    // if(this->_cols != vec.size())
    // {
    //     throw DimensionMismatchException();
    // }
    
    // SymMatrix<DT> ret_matrix(*this);
    // for(int i = 0; i < this->_getArrSize(); )
    // {
    //     for(int j = 0; j < vec.size(); j++)
    //     {
    //         ret_matrix._edges[i] = this->_edges[i] * vec[j];
    //         i++;
    //     }
    // }
    
    // return ret_matrix;
}


//===========================================================MPI SEND/REC================================================================
#ifdef USE_MPI
/*
 * Sends a SymMatrix<T> using MPI_Send.
 * @pram int destination processor
 * @pram int tag 
 * @pram bool sending sparse form, default value is false.
 */
template <typename DT>
void SymMatrix<DT>::MPI_Send_SymMatrix(int dest, int tag, bool sparse = false)
{
    if(sparse)
    {
        MPI_Send(&_SENDING_SYM_SPARSE_FORM, 1 ,MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        std::vector<SparseElement<DT> > sparse_form;
        for (int i = 0; i < this->_size; i++)
        {
            for (int j = i ; j < this->_size; j++)
            {   
                sparse_form.push_back(SparseElement(i,j,(*this)(i,j));
            }
        }
        MPI_Send(&sparse_form.size(), 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<DT>), MPI_BYTE, dest, tag + 3, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&_SENDING_SYM_DENSE_FORM, 1 ,MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        MPI_Send(&this->_size, 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(this->_edges, this->_getArrSize()*sizeof(DT), MPI_BYTE, dest, tag + 3, MPI_COMM_WORLD);
    }
}

/*
 * Sends a SymMatrix<T> using MPI_Bcast.
 * @pram int sourse sender's ID
 * @pram bool sending sparse form, default value is false.
 */
template <typename DT>
void SymMatrix<DT>::MPI_Bcast_Send_SymMatrix(int source, bool sparse = false)
{
    if(sparse)
    {
        MPI_Bcast(&_SENDING_SYM_SPARSE_FORM, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        std::vector<SparseElement<DT> > sparse_form;
        for (int i = 0; i < this->_size; i++)
        {
            for (int j = i ; j < this->_size; j++)
            {   
                sparse_form.push_back(SparseElement(i,j,(*this)(i,j));
            }
        }
        MPI_Bcast(&sparse_form.size(), 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<DT>), MPI_BYTE, source, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&_SENDING_SYM_DENSE_FORM, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_size, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(this->_edges, this->_getArrSize()*sizeof(DT), MPI_BYTE, source, MPI_COMM_WORLD);
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

/*
 * Overloaded = operator copies the content of another matrix to this
 * @pram: SymMatrix<T> 
 */
template <typename DT>
inline SymMatrix<DT>& SymMatrix<DT>::operator=(const SymMatrix<DT>& matrix)
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

    for(int i = 0; i < _getArrSize(); i++)
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
inline DenseMatrix1D<T> DenseMatrix1D<T>::operator+(const DenseMatrix1D<T>& other_matrix) const
{
    // DenseMatrix1D<T> ret_matrix(this->_rows,this->_cols);
    // for(int i = 0; i < this->_getArrSize() ;i++)
    // {   
    //     ret_matrix._edges[i] = this->_edges[i] + other_matrix._edges[i];
    // }
    // return ret_matrix;
}

/*
 * Returns a DenseMatrix1D that is the difference of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::operator-(const DenseMatrix1D<T>& other_matrix) const
{
    // DenseMatrix1D<T> ret_matrix(this->_rows,this->_cols);
    // for(int i = 0; i < this->_getArrSize() ;i++)
    // {   
    //     ret_matrix._edges[i] = this->_edges[i] - other_matrix._edges[i];
    // }
    // return ret_matrix;
}

/*
 * Returns a DenseMatrix1D that is product of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::operator*(const DenseMatrix1D<T>& other_matrix) const
{
    // DenseMatrix2D<T> ret_matrix(this->_rows,this->_cols);
    // T ret_val;
    // for(int i = 0; i < this->_rows;i++)
    // {
    //     for(int j = 0; j < other_matrix._cols;j++)
    //     {
    //         ret_val = 0;
    //         for(int k = 0; k < this->_cols;k++)
    //         {
    //            ret_val += (*this)(i, k) * other_matrix(k, j);
    //         }
    //         ret_matrix(i, j) = ret_val;
    //     }
    // }
    // return ret_matrix;
}

/*
 * Overloaded ostream operator for printing a matrix
 * @pram: std::ostream 
 * @pram: SymMatrix<T>
 */
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

//===========================================================PRIVATE=================================================================
/*
 * Make a deep copy of a matrix object
 * @pram: SymMatrix<T> 
 */
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

/*
 * makes a 1D array of size _rows*_cols
 * @pram: bool fill: if true: initialize values to 0.
 */  
template <typename DT>
inline void SymMatrix<DT>::_initializeMatrix()
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
template<typename DT>
inline size_t SymMatrix<DT>::_getArrSize() const
{
    return (int)(0.5 * this->_size * (this->_size + 1));
}

//===================================================================================================================================

#endif
