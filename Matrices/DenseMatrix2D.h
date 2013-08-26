/*********************************************************************************
 * Dense Matrix Data Structure. This structure uses a two dimensional array that *
 * is dynamically allocated to hold all the values in a matrix.                  * 
 *                                                                               *
 *********************************************************************************/

#ifndef _DenseMatrix_h
#define _DenseMatrix_h

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "MatrixExceptions.h"
#include "SparseElement.h"

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
class DenseMatrix2D;

template <typename T>
std::ostream& operator<< (std::ostream&, const DenseMatrix2D<T>&);

/*
 * DenseMatrix2D class definition and method declarations.
 */
template <typename T>
class DenseMatrix2D
{
private:
    /*
     * Private methods (used internally)
     */
    void _copy(const DenseMatrix2D<T>&);
    void _initializeMatrix(bool);

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
    unsigned int _rows;
    unsigned int _cols;
    T** _edges;
    
public:
    /**************
     *Constructors*
     **************/
    DenseMatrix2D(bool fill = true);
    DenseMatrix2D(int rows, int cols, bool fill = true);
    DenseMatrix2D(const std::string &file_path);
    DenseMatrix2D(const DenseMatrix2D<T>&);

    #ifdef USE_MPI
    DenseMatrix2D(int,MPI_Status&,bool);
    DenseMatrix2D(int,int,MPI_Status&, bool);
    #endif

    /************
     *Destructor*
     ************/
    virtual ~DenseMatrix2D();

    /***********
     *ACCESSORS*
     ***********/
    bool isSquare() const;
    bool isSymmetric() const;
    int getNumberOfRows() const;
    int getNumberOfColumns() const;
    std::vector<SparseElement<T> >getSparseForm() const;
    DenseMatrix2D<T> getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int>& vec_B) const;

    /**********
     *MUTATORS*
     **********/
  
    /**********
    *OPERATIONS*
    **********/
    T getFrobNorm() const;
    T* getTopEigenVector() const;
    DenseMatrix2D<T> transpose() const;
    std::vector<T> getSumOfRows() const;
    std::vector<int> getNeighbors(int vertex) const;
    DenseMatrix2D<T> kron(const DenseMatrix2D<T>& matrix) const;
    DenseMatrix2D<T> diagonalVectorTimesMatrix(const std::vector<T>&) const;
    DenseMatrix2D<T> matrixTimesDiagonalVector(const std::vector<T>&) const;

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
    void operator= (const DenseMatrix2D<T>&);
    bool operator==(const DenseMatrix2D<T>&); 
    DenseMatrix2D<T> operator+(const DenseMatrix2D<T>& other_matrix) const;
    DenseMatrix2D<T> operator-(const DenseMatrix2D<T>& other_matrix) const;
    DenseMatrix2D<T> operator*(const DenseMatrix2D<T>& other_matrix) const;
    friend std::ostream& operator<< <> (std::ostream& stream, const DenseMatrix2D<T>& matrix);
};

//==========================================================CONSTANTS============================================================
template <typename T>
const int DenseMatrix2D<T>::_DEFAULT_MATRIX_SIZE = 1;
template <typename T>
const T DenseMatrix2D<T>::_DEFAULT_MATRIX_ENTRY = 1;
#ifdef USE_MPI
template <typename T>
const int DenseMatrix2D<T>::_SENDING_DENSE_FORM = 0;
template <typename T>
const int DenseMatrix2D<T>::_SENDING_SPARSE_FORM = 1;
template <typename T>
const int DenseMatrix2D<T>::_SENDING_SYM_DENSE_FORM = 2;
template <typename T>
const int DenseMatrix2D<T>::_SENDING_SYM_SPARSE_FORM = 3;
#endif
//==========================================================CONSTRUCTORS============================================================

/*
 * Default constructor:
 * Construct a matrix of size _DEFAULT_MATRIX_SIZE * _DEFAULT_MATRIX_SIZE initialized to 0.
 * @pram bool fill: fills the matrix with 0's. default value is true
 */
template <typename T>
inline DenseMatrix2D<T>::DenseMatrix2D(bool fill = true)
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
inline DenseMatrix2D<T>::DenseMatrix2D(const std::string &file_path)
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
    
    file_reader >> this->_rows;
    file_reader >> this->_cols;
    file_reader >> tmp_x;      //skip the number of lines entry
    _initializeMatrix(true);
  
    while (!file_reader.eof())
    {
        file_reader >> tmp_x;
        file_reader >> tmp_y;
        this->_edges[tmp_x - 1][tmp_y - 1] = _DEFAULT_MATRIX_ENTRY;
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
inline DenseMatrix2D<T>::DenseMatrix2D(int rows, int cols, bool fill = true)
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
inline DenseMatrix2D<T>::DenseMatrix2D(int source, int tag, MPI_Status& stat)
{
    if (isSymmetric)
    {
        MPI_Recv(&this->_rows, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
        this->_cols = _rows;
        _initializeMatrix(false);
        int recv_edges_size = (int)(0.5 * this->_rows * (this->_rows + 1));
        T* recv_edges = new T[recv_edges_size];
        MPI_Recv(recv_edges, recv_edges_size*sizeof(T), MPI_BYTE, source, tag + 2, MPI_COMM_WORLD, &stat);

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
    else
    {
        MPI_Recv(&this->_rows, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
        MPI_Recv(&this->_cols, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
        _initializeMatrix(false);
        int recv_edges_size = this->_rows * this->_cols;
        T* recv_edges = new T[recv_edges_size];
        MPI_Recv(recv_edges, recv_edges_size*sizeof(T), MPI_BYTE, source, tag + 2, MPI_COMM_WORLD, &stat);
        for (int i = 0; i < recv_edges_size; i++)
        {
            this->_edges[i/this->_cols][i%this->_cols] = recv_edges[i];
        }
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
inline DenseMatrix2D<T>::DenseMatrix2D(int source, MPI_Status& stat)
{
    if (isSymmetric)
    {
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        this->_cols = _rows;
        _initializeMatrix(false);
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
    else
    {
        MPI_Bcast(&this->_rows, 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Bcast(&this->_cols, 1, MPI_INT, source, MPI_COMM_WORLD);
        _initializeMatrix(false);
        int recv_edges_size = this->_rows * this->_cols;
        T* recv_edges = new T[recv_edges_size];
        MPI_Bcast(recv_edges, recv_edges_size*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD);
        for (int i = 0; i < recv_edges_size; i++)
        {
            this->_edges[i/this->_cols][i%this->_cols] = recv_edges[i];
        }
    }
}
#endif

/*
 * DensMatrix copy constructor:
 * Construct a copy of a matrix.
 * @pram DenseMatrix2D<T>
 */
template <typename T>
inline DenseMatrix2D<T>::DenseMatrix2D(const DenseMatrix2D<T>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
/*
 * DenseMatrix destructor:
 * deleted the arrays used in the matrix. NOTE: user is responsible for deleting the objects in the matrix if they are dynamically allocated.
 */
template <typename T>
inline DenseMatrix2D<T>::~DenseMatrix2D()
{
    if (_edges != NULL)
    {
        for(int i=0; i < this->_rows; i++)
        {
            delete [] _edges[i];
        }
        delete [] _edges;
    }
}

//===========================================================ACCESSORS===============================================================
/*
 * Returns the number of the rows.
 */
template <typename T>
inline int DenseMatrix2D<T>::getNumberOfRows() const
{
    return _rows;
}

/*
 * Returns the number of the columns.
 */
template <typename T>
inline int DenseMatrix2D<T>::getNumberOfColumns() const
{
    return _cols;
}

/*
 * Returns true if the matrix is square (i.e rows == cols)
 */
template <typename T>
inline bool DenseMatrix2D<T>::isSquare() const
{
    return (this->_rows == this->_cols);
}

/*
 * Returns true if the matrix is symmetric (i.e values at i, j are equal to values at j, i)
 */
template <typename T>
inline bool DenseMatrix2D<T>::isSymmetric() const
{
    //sym. matrix has to be square
    if (this->_rows != this->_cols)
    {
        return false;
    }

    //checking for entries to be equal
    for(int i=0; i < this->_rows; i++)
    {
        for(int j=0; j<this->_cols;j++)
        {
            if (this->_edges[i][j] != this->_edges[j][i])
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
inline std::vector<SparseElement<T> > DenseMatrix2D<T>::getSparseForm() const
{    
    std::vector<SparseElement<T> > sparse_form;
    for(int i = 0; i < this->_rows; i++)
    {
        for(int j = 0; j < this->_cols; j++)
        {
            if(this->_edges[i][j] != 0)
            {
                sparse_form.push_back(SparseElement<T>(i,j, this->_edges[i][j]));
            }
        }
    }
    return sparse_form;
}


/*
 * Returns a DenseMatrix2D that has the rows that are marked 1 in vec_A, columns that are marked 1 in vec_B
 * @pram std::vector<int>: vector of 0's and 1's for row selection
 * @pram std::vector<int>: vector of 0's and 1's for column selection
 */
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::getScatteredSelection(const std::vector<int>& vec_A, const std::vector<int>& vec_B) const
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
    DenseMatrix2D<T> res_matrix(num_in_A, num_in_B);
    int counter = 0;
    for (int i=0; i< vec_A.size(); i++)
    {
        for(int j=0; j< vec_B.size(); j++)
        {
            if ( vec_A[i] == 1 && vec_B[j] ==1)
            {
                res_matrix._edges[counter/num_in_B][counter%num_in_B] = this->_edges[i][j];
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
inline T DenseMatrix2D<T>::getFrobNorm() const
{
    T ret_val = 0;

    for(int i = 0; i < this->_rows; i++)
    {
        for(int j = 0; j < this->_cols; j++)
        {
            ret_val+= this->_edges[i][j] * this->_edges[i][j];
        }
    }
    return ret_val;
}

/*
 * Returns a std::vector<T> that contains the values of the eigenvector associated to the largest eigenvalue
 */
template <typename T>
inline T* DenseMatrix2D<T>::getTopEigenVector()
{
#ifdef ARPACK
    int arr_size = (this->_rows*(this->_rows+1))/2;
    T nzval[arr_size];
    int counter=0;
    for(int i=0; i<this->_rows; i++)
    {
        for(int j=i;j < this->_rows ;j++)
        {
            nzval[counter++]= this->_edges[i][j];
        }
    }
    
    ARdsSymMatrix<T> ARMatrix(this->_rows,nzval,'L');
    ARluSymStdEig<T> eigProb(1, ARMatrix, "LM",10);
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
        for (int i=0; i< this->_rows; i++) 
        {
            for(int j = 0; j < this->_cols; j++)
            {
                A_eigen(i,j) = this->_edges[i][j];
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
 * Returns a DenseMatrix2D<T> of the transpose of this.
 */
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::transpose() const
{
    DenseMatrix2D<T> ret_matrix(this->_rows,this->_cols, false);

    for(int i=0;i<this->_rows;i++)
    {
        for(int j=0;j<this->_cols;j++)
        {
            ret_matrix._edges[j][i] = this->_edges[i][j];
        }
    }
    return ret_matrix;
}

/*
 * Returns a std::vector<T> that contains sum of the values in each row
 */
template <typename T>
inline std::vector<T> DenseMatrix2D<T>::getSumOfRows() const
{
    std::vector<T> sum_vec(this->_rows);
    for(int i=0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++)
        {
            sum_vec[i] += this->_edges[i][j];
        }
    }
    return sum_vec;
}

/*
 * Returns a std::vector<int> of the elements with value of 1 in a columns.
 * @pram int vertex
 */
template <typename T>
inline std::vector<int> DenseMatrix2D<T>::getNeighbors(int vertex) const
{
  std::vector<int> neighbors;
  for(int i = 0; i < this->_rows; i++)
  {
    if(this->_edges[i][vertex] == 1)
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
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::kron(const DenseMatrix2D<T>& matrix) const
{
    // also checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    //Initializing and allocating the product matrix
    int prod_size = this->_rows * matrix._rows;
    DenseMatrix2D<T> prod_matrix(prod_size, prod_size, false);
    
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
                    prod_matrix._edges[(i_outer*matrix._rows) + i_inner][(j_outer*matrix._cols)+j_inner] =
                    this->_edges[i_outer][j_outer] * matrix._edges[i_inner][j_inner];
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
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::diagonalVectorTimesMatrix(const std::vector<T>& vec) const
{    
    if(_cols != vec.size())
    {
        throw DimensionMismatchException();
    }
    
    DenseMatrix2D<T> ret_matrix(*this);
    for(int i = 0; i < _rows; i++)
    {
        for(int j = 0; j < _cols; j++ )
        {
            ret_matrix._edges[i][j] = this->_edges[i][j] * vec[i];
        }
    }
    return ret_matrix;
}

/*
 * Returns a DenseMatrix2D<T> of the product of this matrix a vector that contains the diagonal values of a diagonal matrix.
 * @pram std::vector<T> diagonal entires of a diagonal matrix
 */
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::matrixTimesDiagonalVector(const std::vector<T>& vec) const
{    
    if(_cols != vec.size())
    {
        throw DimensionMismatchException();
    }

    DenseMatrix2D<T> ret_matrix(*this);
    for(int i=0; i < _cols; i++)
    {
        for(int j=0; j < _rows; j++)
        {
            ret_matrix._edges[j][i] = this->_edges[j][i]*vec[i];
        }
    }
    return ret_matrix;
}

//===========================================================MPI SEND/REC================================================================
#ifdef USE_MPI
/*
 * Sends a DenseMatrix2D<T> using MPI_Send.
 * @pram int destination processor
 * @pram int tag 
 * @pram bool sending sparse form, default value is false.
 */
template <typename T>
inline void DenseMatrix2D<T>::MPI_Send_Matrix(int dest, int tag, bool sparse = false)
{
    if(sparse)
    {
        MPI_Send(&_SENDING_SPARSE_FORM, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        std::vector<T> sparse_form = this->getSparseForm();
        MPI_Send(&sparse_form.size(), 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
        MPI_Send(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<T>), MPI_BYTE, dest, tag + 3, MPI_COMM_WORLD);

    }
    else
    {
        MPI_Send(&_SENDING_DENSE_FORM, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
        //TODO
    }
}

/*
 * Sends a DenseMatrix2D<T> using MPI_Bcast.
 * @pram int sourse sender's ID
 * @pram bool sending sparse form, default value is false.
 */
template <typename T>
inline void DenseMatrix2D<T>::MPI_Bcast_Send_Matrix(int source, bool sparse = false)
{
    if(sparse)
    {
        MPI_Bcast(&_SENDING_SPARSE_FORM, 1, MPI_INT, source, MPI_COMM_WORLD);
        std::vector<T> sparse_form = this->getSparseForm();
        MPI_Send(&sparse_form.size(), 1, MPI_INT, source, MPI_COMM_WORLD);
        MPI_Send(&sparse_form[0], sparse_form.size()*sizeof(SparseElement<T>), MPI_BYTE, sourse, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&_SENDING_DENSE_FORM, 1, MPI_INT, source, MPI_COMM_WORLD);
        //TODO

    }
}
//==========================================================OPERATORS================================================================

/*
 * Overloaded () operator for accessing and changing the values inside a matrix
 * @pram: int i
 * @pram: int j
 */
template <typename T>
inline T&  DenseMatrix2D<T>::operator()(int i, int j)
{
    return this->_edges[i][j];
}

/*
 * Overloaded = operator copies the content of another matrix to this
 * @pram: DenseMatrix2D<T> 
 */
template <typename T>
inline void DenseMatrix2D<T>::operator=(const DenseMatrix2D<T>& matrix)
{
    if (_edges != NULL)
    {
        for(int i=0; i < this->_rows; i++)
        {
            delete [] _edges[i];
        }
        delete [] _edges;
    }   
    _copy(matrix);
}

/*
 * Overloaded == operator to compare two matrices
 * @pram: DenseMatrix2D<T> 
 */
template <typename T>
inline bool DenseMatrix2D<T>::operator==(const DenseMatrix2D<T>& matrix)
{
    if ((this->_cols != matrix._cols) || (this->_rows != matrix._rows))
    {
        return false;
    }

    for (int i = 0; i < this->_rows; i++)
    {
        for (int j = 0; i < this->_cols; j++)
        {
            if ( this->_edges[i][j] != matrix._edges[i][j])
            {
                return false;
            }
        }   
    }
    return true;
}

/*
 * Returns a DenseMatrix2D that is the sum of this and other_matrix
 * @pram: DenseMatrix2D<T> 
 */
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::operator+(const DenseMatrix2D<T>& other_matrix) const
{
  DenseMatrix2D<T> ret_matrix(this->_rows,this->_cols, false);
  for(int i = 0; i < this->_rows ;i++)
  {   
      for(int j = 0; j < this->_cols; j++)
      {
           ret_matrix._edges[i][j] = this->_edges[i][j] + other_matrix._edges[i][j];
      }
  }
  return ret_matrix;
}

/*
 * Returns a DenseMatrix2D that is the difference of this and other_matrix
 * @pram: DenseMatrix2D<T> 
 */
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::operator-(const DenseMatrix2D<T>& other_matrix) const
{
  DenseMatrix2D<T> ret_matrix(this->_rows,this->_cols, false);
  for(int i = 0; i < this->_rows ;i++)
  {   
      for(int j = 0; j < this->_cols; j++)
      {
           ret_matrix._edges[i][j] = this->_edges[i][j] - other_matrix._edges[i][j];
      }
  }
  return ret_matrix;
}

/*
 * Returns a DenseMatrix2D that is product of this and other_matrix
 * @pram: DenseMatrix2D<T> 
 */
template <typename T>
inline DenseMatrix2D<T> DenseMatrix2D<T>::operator*(const DenseMatrix2D<T>& other_matrix) const
{
    DenseMatrix2D<T> ret_matrix(this->_rows,this->_cols, false);
    T ret_val;
    for(int i = 0; i < this->_rows;i++)
    {
        for(int j = 0; j < other_matrix._cols;j++)
        {
            ret_val = 0;
            for(int k = 0; k < this->_cols;k++)
            {
               ret_val += this->_edges[i][k] * other_matrix._edges[k][j];
            }
            ret_matrix._edges[i][j] = ret_val;
        }
    }
  
    return ret_matrix;
}

/*
 * overloaded ostream operator for printing a matrix
 * @pram: std::ostream 
 * @pram: DenseMatrix2D<T>
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const DenseMatrix2D<T>& matrix)
{
    stream<< "Size: " << matrix._rows << "*" << matrix._cols << '\n';
    for (int i = 0; i < matrix._rows; i++)
    {
        for( int j = 0; j< matrix._cols; j++)
        {
            stream << matrix._edges[i][j] << ' ';
        }
        stream << "\n";
    }
    stream << "\n\n\n";
    return stream;
}

//===========================================================PRIVATE=================================================================
/*
 * Make a deep copy of a matrix object
 * @pram: DenseMatrix2D<T> 
 */
template <typename T>
inline void DenseMatrix2D<T>::_copy(const DenseMatrix2D<T>& matrix)
{
    this->_rows = matrix._rows;
    this->_cols = matrix._cols;
    _initializeMatrix(false);

    for(int i=0; i < this->_rows; i++)
    {
        memcpy(this->_edges[i], matrix._edges[i], this->_cols * sizeof(T));
    }
}

/*
 * makes a 2D array of size _rows*_cols
 * @pram: bool fill: if true: initialize values to 0.
 */
template <typename T>
inline void DenseMatrix2D<T>::_initializeMatrix(bool fill)
{
    try
    {
        this->_edges = new T*[this->_rows];
      
        if (fill)
        {
            for(int i=0; i < this->_rows; i++)
            {
                this->_edges[i] = new T[this->_cols]();
            }
        }
        else
        {
            for(int i=0; i < this->_rows; i++)
            {
                this->_edges[i] = new T[this->_cols];
            }
        }
    }
    catch (std::bad_alloc& e)
    {
        throw OutOfMemoryException();
    }
}

//===================================================================================================================================
#endif