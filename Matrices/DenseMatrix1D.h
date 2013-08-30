/*********************************************************************************
 * Dense Matrix Data Structure. This structure uses a one dimensional array that *
 * is dynamically allocated to hold all the values in a matrix.                  * 
 *                                                                               *
 *********************************************************************************/

#ifndef _DenseMatrix1D_h
#define _DenseMatrix1D_h

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "MatrixExceptions.h"
#include "SparseElement.h"

#ifdef ARPACK
#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"
#endif

#ifdef EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

template <typename T>
class DenseMatrix1D;

template <typename T>
std::ostream& operator<< (std::ostream&, const DenseMatrix1D<T>&);

/*
 * DenseMatrix1D class definition and method declarations.
 */
template <typename T>
class DenseMatrix1D
{
private:
    /*
     * Private methods (used internally)
     */
    size_t _getArrSize() const;
    void _initializeMatrix(bool);
    void _copy(const DenseMatrix1D<T>&);

    /*
     * Private constants.
     */
    static const int _DEFAULT_MATRIX_SIZE;
    static const  T _DEFAULT_MATRIX_ENTRY;
#ifdef USE_MPI
    static const int _SENDING_SPARSE_FORM;
    static const int _SENDING_DENSE_FORM;
    static const int _SENDING_SYM_SPARSE_FORM;
    static const int _SENDING_SYM_DENSE_FORM;
#endif
protected:
    int _rows;
    int _cols;
    T* _edges;
    
    
public:
    /**************
     *Constructors*
     **************/
    DenseMatrix1D(bool fill = true);
    DenseMatrix1D(int, int, bool fill = true);
    DenseMatrix1D(const DenseMatrix1D<T>&);
    DenseMatrix1D(const std::string&);

    #ifdef USE_MPI
    DenseMatrix1D(int,MPI_Status&);
    DenseMatrix1D(int,int,MPI_Status&);
    #endif
    
    /************
     *Destructor*
     ************/
    virtual ~DenseMatrix1D();
    
    /***********
     *ACCESSORS*
     ***********/
    bool isSquare();
    bool isSymmetric();
    int getNumberOfRows();
    int getNumberOfColumns();
    std::vector<SparseElement<T> >getSparseForm();
    DenseMatrix1D<T> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B);
    
    /**********
     *MUTATORS*
     **********/

    /**********
    *OPERATIONS*
    **********/
    T getFrobNorm();
    T* getTopEigenVector();
    DenseMatrix1D<T> transpose();
    std::vector<T> getSumOfRows();
    std::vector<int> getNeighbors(int vertex);
    DenseMatrix1D<T> kron(const DenseMatrix1D<T>& matrix);
    DenseMatrix1D<T> diagonalVectorTimesMatrix(const std::vector<T>&);
    DenseMatrix1D<T> matrixTimesDiagonalVector(const std::vector<T>&);

    /*************
    *  MPI Send  *
    **************/
    #ifdef USE_MPI
    void MPI_Send_Matrix(int, int,bool sparse = false);
    void MPI_Bcast_Send_Matrix(int,bool sparse = false);
    #endif

    /**********
     *OPERATORS*
     **********/
    T& operator()(int i, int j);
    void operator= (const DenseMatrix1D<T>&);
    bool operator==(const DenseMatrix1D<T>&);
    DenseMatrix1D<T> operator+(const DenseMatrix1D<T>& other_matrix);
    DenseMatrix1D<T> operator-(const DenseMatrix1D<T>& other_matrix);
    DenseMatrix1D<T> operator*(const DenseMatrix1D<T>& other_matrix);
    friend std::ostream& operator<< <> (std::ostream& stream, const DenseMatrix1D<T>& matrix);
};

//==========================================================CONSTANTS============================================================
template <typename T>
const int DenseMatrix1D<T>::_DEFAULT_MATRIX_SIZE = 1;
template <typename T>
const T DenseMatrix1D<T>::_DEFAULT_MATRIX_ENTRY = 1;
#ifdef USE_MPI
template <typename T>
const int DenseMatrix1D<T>::_DENSE_FORM = 0;
template <typename T>
const int DenseMatrix1D<T>::_SPARSE_FORM = 1;
template <typename T>
const int DenseMatrix1D<T>::_SYM_DENSE_FORM = 2;
template <typename T>
const int DenseMatrix1D<T>::_SYM_SPARSE_FORM = 3;
#endif
//==========================================================CONSTRUCTORS============================================================
/*
 * Default constructor:
 * Construct a matrix of size _DEFAULT_MATRIX_SIZE * _DEFAULT_MATRIX_SIZE initialized to 0.
 * @pram bool fill: fills the matrix with 0's. default value is true
 */
template <typename T>
inline DenseMatrix1D<T>::DenseMatrix1D(bool fill)
{
    this->_rows = _DEFAULT_MATRIX_SIZE;
    this->_cols = _DEFAULT_MATRIX_SIZE;
    _initializeMatrix(fill);
}

/*
 * DensMatrix constructor:
 * Construct a matrix by reading a matrix file, specified in readme.txt file the nodes that exist will have _DEFAULT_MATRIX_ENTRY value.
 * @pram std::string : path to the file
 */
template<typename T>
inline DenseMatrix1D<T>::DenseMatrix1D(const std::string& file_path)
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
    
    _initializeMatrix(true);
    file_reader >> tmp_x;      //skip the line number
    
    while (!file_reader.eof())
    {
        file_reader >> tmp_x;
        file_reader >> tmp_y;
        (*this)(tmp_x - 1 ,tmp_y - 1) = _DEFAULT_MATRIX_ENTRY;
    }
    
    file_reader.close();
}

/*
 * DensMatrix constructor:
 * Construct a matrix n*m matrix initialized to 0.
 * @pram int rows: number of rows
 * @pram int cols: number of columns
 * @pram bool fill: fills the matrix with 0's. default value is true
 */
template <typename T>
inline DenseMatrix1D<T>::DenseMatrix1D(int rows, int cols, bool fill)
{
    this->_rows = rows;
    this->_cols = cols;
    _initializeMatrix(fill);
}

#ifdef USE_MPI
/*
 * DensMatrix constructor:
 * Construct a matrix by receiving matrix data from a MPI_Send_Matrix call.
 * Note the matrix being send should be of type DenseMatrix1D, DenseMatrix2D, SymmMatrix
 * @pram int source: Sender's ID
 * @pram int tag: sender's tag
 * @pram MPI_Status: MPI_Status object
 */
template <typename T>
inline DenseMatrix1D<T>::DenseMatrix1D(int source, int tag, MPI_Status& stat)
{
    int recv_format;
    MPI_Recv(&recv_format, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
    if(recv_format == _SYM_DENSE_FORM)
    {
        MPI_Recv(&this->_rows, 1, MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
        this->_cols = _rows;
        _initializeMatrix(false);
        int recv_edges_size = (int)(0.5 * this->_rows * (this->_rows + 1));
        T* recv_edges = new T[recv_edges_size];
        MPI_Recv(recv_edges, recv_edges_size*sizeof(T), MPI_BYTE, source, tag + 3, MPI_COMM_WORLD, &stat);

        int counter = 0;
        for (int i = 0; i < this->_rows; i++)
        {
            for (int j = 0; j <= i; j++)
            {
               this->_edges[j][i] = this->_edges[i][j] = recv_edges[counter++];  
            }
        }
        delete [] recv_edges;
    }
    else if (recv_format == _DENSE_FORM)
    {
        MPI_Recv(&this->_rows, 1, MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
        MPI_Recv(&this->_cols, 1, MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
        _initializeMatrix(false);
        MPI_Recv(this->_edges, _getArrSize()*sizeof(T), MPI_BYTE, source, tag + 4, MPI_COMM_WORLD, &stat);
    }
    else if (recv_format == _SYM_SPARSE_FORM)
    {
        MPI_Recv(&this->_rows, 1, MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
        this->_cols = _rows;
        _initializeMatrix(true);
        int recv_edges_size;
        MPI_Recv(&recv_edges_size, 1, MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
        SparseElement<T>* recv_edges = new SparseElement<T>[recv_edges_size];
        MPI_Recv(recv_edges, recv_edges_size*sizeof(SparseElement<T>), MPI_BYTE, source, tag + 4, MPI_COMM_WORLD, &stat);
        for(int i = 0;  i < recv_edges_size; i++)
        {
            (*this)(recv_edges[i].getJ(),recv_edges[i].getI()) = (*this)(recv_edges[i].getI(), recv_edges[i].getJ()) = recv_edges[i].getValue();
        }
        delete recv_edges;
    }
    else if (recv_format == _SPARSE_FORM)
    {
        MPI_Recv(&this->_rows, 1, MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
        MPI_Recv(&this->_cols, 1, MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
        _initializeMatrix(true);
        int recv_edges_size;
        MPI_Recv(&recv_edges_size, 1, MPI_INT, source, tag + 4, MPI_COMM_WORLD, &stat);
        SparseElement<T>* recv_edges = new SparseElement<T>[recv_edges_size];
        MPI_Recv(recv_edges, recv_edges_size*sizeof(SparseElement<T>), MPI_BYTE, source, tag + 5, MPI_COMM_WORLD, &stat);
        for(int i = 0;  i < recv_edges_size; i++)
        {
            (*this)(recv_edges[i].getI(), recv_edges[i].getJ()) = recv_edges[i].getValue();
        }
        delete recv_edges;
    }
}

/*
 * DensMatrix constructor:
 * Construct a matrix by receiving matrix data from a MPI_Bcast_Send_Matrix call.
 * Note the matrix being send should be of type DenseMatrix1D, DenseMatrix2D, SymmMatrix
 * @pram int source: Sender's ID
 * @pram MPI_Status: MPI_Status object
 */
template <typename T>
inline DenseMatrix1D<T>::DenseMatrix1D(int source, MPI_Status& stat)
{
    int recv_format;
    MPI_Bcast(&recv_format, 1, MPI_INT, source, MPI_COMM_WORLD, &stat);
    if(recv_format == _SYM_DENSE_FORM)
    {
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        this->_cols = _rows;
        _initializeMatrix(true);
        int recv_edges_size = (int)(0.5 * this->_rows * (this->_rows + 1));
        T* recv_edges = new T[recv_edges_size];
        MPI_Bcast(recv_edges, recv_edges_size*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD);

        int counter = 0;
        for (int i = 0; i < this->_rows; i++)
        {
            for (int j = 0; j <= i; j++)
            {
               this->_edges[j][i] = this->_edges[i][j] = recv_edges[counter++];  
            }
        }
        delete [] recv_edges;
    }
    else if (recv_format == _DENSE_FORM)
    {
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_cols, 1, MPI_INT, source, MPI_COMM_WORLD);
        _initializeMatrix(false);
        MPI_Bcast(this->_edges, _getArrSize()*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD, &stat);
    }
    else if (recv_format == _SYM_SPARSE_FORM)
    {
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        this->_cols = _rows;
        _initializeMatrix(true);
        int recv_edges_size;
        MPI_Bcast(&recv_edges_size, 1, MPI_INT, source, MPI_COMM_WORLD);
        T* recv_edges = new T[recv_edges_size];
        MPI_Bcast(recv_edges, recv_edges_size*sizeof(SparseElement<T>), MPI_BYTE, source, MPI_COMM_WORLD);
        for(int i = 0;  i < recv_edges_size; i++)
        {
            (*this)(recv_edges[i].getJ(),recv_edges[i].getI()) = (*this)(recv_edges[i].getI(), recv_edges[i].getJ()) = recv_edges[i].getValue();
        }
    }
    else if (recv_format == _SPARSE_FORM)
    {
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        this->_cols = _rows;
        _initializeMatrix(true);
        int recv_edges_size;
        MPI_Bcast(&recv_edges_size, 1, MPI_INT, source, MPI_COMM_WORLD);
        T* recv_edges = new T[recv_edges_size];
        MPI_Bcast(recv_edges, recv_edges_size*sizeof(SparseElement<T>), MPI_BYTE, source, MPI_COMM_WORLD);
        for(int i = 0;  i < recv_edges_size; i++)
        {
            (*this)(recv_edges[i].getI(), recv_edges[i].getJ()) = recv_edges[i].getValue();
        }

    }
}
#endif

/*
 * DensMatrix copy constructor:
 * Construct a copy of a matrix.
 * @pram DenseMatrix1D<T>
 */
template <typename T>
inline DenseMatrix1D<T>::DenseMatrix1D(const DenseMatrix1D<T>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
/*
 * DenseMatrix destructor:
 * deleted the arrays used in the matrix. NOTE: user is responsible for deleting the objects in the matrix if they are dynamically allocated.
 */
template <typename T>
inline DenseMatrix1D<T>::~DenseMatrix1D()
{
    if (this->_edges != NULL)
    {
        delete [] this->_edges;
    }
}

//===========================================================ACCESSORS===============================================================
/*
 * Returns the number of the rows.
 */
template <typename T>
inline int DenseMatrix1D<T>::getNumberOfRows()
{
    return this->_rows;
}

/*
 * Returns the number of the columns.
 */
template <typename T>
inline int DenseMatrix1D<T>::getNumberOfColumns()
{
    return this->_cols;
}

/*
 * Returns true if the matrix is square (i.e rows == cols)
 */
template <typename T>
inline bool DenseMatrix1D<T>::isSquare()
{
    return (this->_rows == this->_cols);
}

/*
 * Returns true if the matrix is symmetric (i.e values at i, j are equal to values at j, i)
 */
template <typename T>
inline bool DenseMatrix1D<T>::isSymmetric()
{
    //sym. matrix has to be square
    if (this->_rows != this->_cols)
    {
        return false;
    }
    
    //checking for entries to be equal
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

/*
 * Returns a std::vector<SparseElement<T>> of SparseElement objects that contain the i, j and value of the none-zero edges.
 */
template <typename T>
inline std::vector<SparseElement<T> > DenseMatrix1D<T>::getSparseForm()
{    
    std::vector<SparseElement<T> > sparse_form;
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        if(this->_edges[i] != 0)
        {
            sparse_form.push_back(SparseElement<T>(i/this->_cols, i%this->_cols, this->_edges[i]));
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
inline DenseMatrix1D<T> DenseMatrix1D<T>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int> vec_B)
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
    DenseMatrix1D<T> res_matrix(num_in_A, num_in_B);
    
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



//===========================================================MUTATORS================================================================

//===========================================================OPERATIONS================================================================
/*
 * Returns the frobenius norm.
 */
template <typename T>
inline T DenseMatrix1D<T>::getFrobNorm()
{
    T ret_val = 0;
    for (int i = 0; i < this->_getArrSize(); i++)
    {
        ret_val += this->_edges[i] * this->_edges[i];
    }

    return ret_val;
}

/*
 * Returns a std::vector<T> that contains the values of the eigenvector associated to the largest eigenvalue
 */
template <typename T>
inline T* DenseMatrix1D<T>::getTopEigenVector()
{
#ifdef ARPACK
    int arr_size = (0.5 * this->_rows * (this->_rows+1));
    T sym_edges[arr_size];
    int counter = 0;
    
    for(int i = 0; i < this->_rows; i++)
    {
        for(int j = i; j < this->_cols; j++)
        {
            sym_edges[counter] = this->_edges[i * this->_cols +j];
            counter++;
        }
    }
    
    ARdsSymMatrix<T> ARMatrix(this->_rows, sym_edges, 'L');
    ARluSymStdEig<T> eigProb(1, ARMatrix, "LM", 10);
    eigProb.FindEigenvectors();
    T* eigen_vector = new T[this->_rows];
     
    for (int i=0; i < this->_rows ; i++)
    {
        eigen_vector[i] = eigProb.Eigenvector(0,i);
    }
    return eigen_vector;
#endif

#ifdef EIGEN
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> A_eigen = Eigen::MatrixXd::Zero(this->_rows, this->_cols);
    for (long i=0; i< this->_rows; i++) 
    {
        for(int j = 0; j < this->_cols; j++)
        {
            A_eigen(i,j) = this->_edges[i * this->_cols +j];
        }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(A_eigen);
    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> evals_eigen = es.eigenvalues();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> evecs_eigen = es.eigenvectors();
    T* eigen_vector = new T[this->_rows];

    for ( int i=0; i < this->_cols; i++)
    {   
        if (evals_eigen(i) == evals_eigen.maxCoeff())
        {
            for (int j= 0; j < this->_rows ; j++)
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
 * Returns a DenseMatrix1D<T> of the transpose of this.
 */
template <typename T>
DenseMatrix1D<T> DenseMatrix1D<T>::transpose()
{
    DenseMatrix1D<T> ret_matrix(this->_rows,this->_cols);

    for(int i=0;i<this->_rows;i++)
    {
        for(int j=0;j<this->_cols;j++)
        {
            ret_matrix._edges[j * ret_matrix._cols + i] = this->_edges[i *this->_cols + j];
        }
    }
    return ret_matrix;
}

/*
 * returns a vector of type T where each entry corresponds to sum of entries in a matrix row.
 * Delete this object after usage.
 * @return std::vector<T>
 */
template <typename T>
inline std::vector<T> DenseMatrix1D<T>::getSumOfRows()
{
    std::vector<T> sum_vector(this->_rows);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        sum_vector[(i / this->_rows)] += this->_edges[i];
    }
    
    return sum_vector;
}

/*
 * Returns a std::vector<int> of the elements with value of 1 in a columns.
 * @pram int vertex
 */
template <typename T>
inline std::vector<int> DenseMatrix1D<T>::getNeighbors(int vertex)
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

/*
 * Returns a DenseMatrix1D that is the kronecker product of this and another matrix.
 * @pram DenseMatrix1D
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::kron(const DenseMatrix1D<T>& matrix)
{
    // checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    //Initializing and allocating the product matrix
    int prod_size = this->_rows * matrix._rows;
    DenseMatrix1D<T> prod_matrix(prod_size, prod_size);
    
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
                    prod_matrix._edges[((i_outer*matrix._rows) + i_inner)* prod_matrix._cols + (j_outer*matrix._cols) + j_inner] 
                        = matrix._edges[i_inner*matrix._cols+ j_inner] * this->_edges[i_outer*this->_cols + j_outer];
                }
            }
        }
    }
    return prod_matrix;
}

/*
 * Returns a DenseMatrix1D<T> of a vector that contains the diagonal values of a diagonal matrix times this matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::diagonalVectorTimesMatrix(const std::vector<T>& vec)
{
    if(this->_rows != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    DenseMatrix1D<T> ret_matrix(*this);
    for(int i = 0; i < this->_getArrSize(); i++)
    {
        ret_matrix._edges[i] = vec[i/this->_rows] * this->_edges[i];
    }
    
    return ret_matrix;
}

/*
 * Returns a DenseMatrix1D<T> of the product of this matrix a vector that contains the diagonal values of a diagonal matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::matrixTimesDiagonalVector(const std::vector<T>& vec)
{
    if(this->_cols != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    DenseMatrix1D<T> ret_matrix(*this);
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

//===========================================================MPI SEND/REC================================================================
#ifdef USE_MPI
/*
 * Sends a DenseMatrix1D<T> using MPI_Send.
 * @pram int destination processor
 * @pram int tag 
 * @pram bool sending sparse form, default value is false.
 */
template <typename T>
inline void DenseMatrix1D<T>::MPI_Send_Matrix(int dest, int tag, bool sparse)
{
    if(sparse)
    {
        MPI_Send(&_SPARSE_FORM, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        MPI_Send(&this->_rows, 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(&this->_cols, 1, MPI_INT, dest, tag + 3, MPI_COMM_WORLD);
        std::vector<T> sparse_form = this->getSparseForm();
        MPI_Send(&sparse_form.size(), 1, MPI_INT, dest, tag + 4, MPI_COMM_WORLD);
        MPI_Send(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<T>), MPI_BYTE, dest, tag + 5, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&_DENSE_FORM, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        MPI_Send(&this->_rows, 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(&this->_cols, 1, MPI_INT, dest, tag + 3, MPI_COMM_WORLD);
        MPI_Send(this->_edges, this->_getArrSize()*sizeof(T), MPI_BYTE, dest, tag + 4, MPI_COMM_WORLD);
    }
}

/*
 * Sends a DenseMatrix1D<T> using MPI_Bcast.
 * @pram int sourse sender's ID
 * @pram bool sending sparse form, default value is false.
 */
template <typename T>
inline void DenseMatrix1D<T>::MPI_Bcast_Send_Matrix(int source, bool sparse)
{
    if(sparse)
    {
        MPI_Bcast(&_SPARSE_FORM, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_cols, 1, MPI_INT, source, MPI_COMM_WORLD);
        std::vector<T> sparse_form = this->getSparseForm();
        MPI_Send(&sparse_form.size(), 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Send(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<T>), MPI_BYTE, source, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&_SENDING_DENSE_FORM, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_cols, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_size, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(this->_edges, this->_getArrSize()*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD);
    }
}
#endif
//==========================================================OPERATORS================================================================

/*
 * Overloaded () operator for accessing and changing the values inside a matrix
 * @pram: int i
 * @pram: int j
 */
template <typename T>
inline T& DenseMatrix1D<T>::operator()(int i, int j)
{
    return this->_edges[(i * this->_cols) + j];
}

/*
 * Overloaded = operator copies the content of another matrix to this
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline void DenseMatrix1D<T>::operator=(const DenseMatrix1D<T>& matrix)
{
    if (this->_edges != NULL)
    {
        delete [] this->_edges;
    }
    _copy(matrix);
}

/*
 * Overloaded == operator to compare two matrices
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline bool DenseMatrix1D<T>::operator==(const DenseMatrix1D<T>& matrix)
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

/*
 * Returns a DenseMatrix1D that is the sum of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::operator+(const DenseMatrix1D<T>& other_matrix)
{
    DenseMatrix1D<T> ret_matrix(this->_rows,this->_cols);
    for(int i = 0; i < this->_getArrSize() ;i++)
    {   
        ret_matrix._edges[i] = this->_edges[i] + other_matrix._edges[i];
    }
    return ret_matrix;
}

/*
 * Returns a DenseMatrix1D that is the difference of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::operator-(const DenseMatrix1D<T>& other_matrix)
{
    DenseMatrix1D<T> ret_matrix(this->_rows,this->_cols);
    for(int i = 0; i < this->_getArrSize() ;i++)
    {   
        ret_matrix._edges[i] = this->_edges[i] - other_matrix._edges[i];
    }
    return ret_matrix;
}

/*
 * Returns a DenseMatrix1D that is product of this and other_matrix
 * @pram: DenseMatrix1D<T> 
 */
template <typename T>
inline DenseMatrix1D<T> DenseMatrix1D<T>::operator*(const DenseMatrix1D<T>& other_matrix)
{
    DenseMatrix1D<T> ret_matrix(this->_rows,this->_cols);
    T ret_val;
    for(int i = 0; i < this->_rows;i++)
    {
        for(int j = 0; j < other_matrix._cols;j++)
        {
            ret_val = 0;
            for(int k = 0; k < this->_cols;k++)
            {
               ret_val += this->_edges[i*this->_cols + k] * other_matrix._edges[k*other_matrix._cols + j];
            }
            ret_matrix._edges[i * ret_matrix._cols + j] = ret_val;
        }
    }
    return ret_matrix;
}

/*
 * overloaded ostream operator for printing a matrix
 * @pram: std::ostream 
 * @pram: DenseMatrix1D<T>
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const DenseMatrix1D<T>& matrix)
{
    stream<< "Size: " << matrix._rows << "*" << matrix._cols << '\n';
    for (int i = 0; i < matrix._rows; i++)
    {
        for( int j = 0; j< matrix._cols; j++)
        {
            stream << matrix._edges[i*matrix._cols +j] << ' ';
        }
        stream << "\n";
    }
    stream << "\n\n\n";
    return stream;
}

//===========================================================PRIVATE=================================================================
/*
 * Make a deep copy of a matrix object
 * @pram: DenseMatrix1<T> 
 */
template <typename T>
inline void DenseMatrix1D<T>::_copy(const DenseMatrix1D<T>& matrix)
{
    this->_rows = matrix._rows;
    this->_cols = matrix._cols;
    _initializeMatrix(false);
    memcpy(this->_edges, matrix._edges, this->_getArrSize() * sizeof(T));
}

/*
 * makes a 1D array of size _rows*_cols
 * @pram: bool fill: if true: initialize values to 0.
 */  
template <typename T>
inline void DenseMatrix1D<T>::_initializeMatrix(bool fill)
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
inline size_t DenseMatrix1D<T>::_getArrSize() const
{
    return (this->_rows * this->_cols);
}
    
//===================================================================================================================================
    
    
#endif
    