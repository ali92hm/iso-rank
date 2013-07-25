//
//  DenseMatrix.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _DenseMatrix_h
#define _DenseMatrix_h

#define DEFAULT_MATRIX_SIZE 1

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
#include <stdexcept>
#include <limits>
#include <math.h>
#include "MatrixExceptions.h"
#include "../util.h"

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


template <typename DT>
class DenseMatrix;

template <typename DT>
struct sparse_matrix_element{
    int row_index;
    int col_index;
    DT value;
};


template <typename DT>
std::ostream& operator<< (std::ostream&, const DenseMatrix<DT>&);

template <typename DT>
class DenseMatrix
{
private:
    std::ifstream _file_reader;
    void _copy(const DenseMatrix<DT>&);
    bool _initilalizeMatrix();
    static void _split(const std::string &line, std::vector<unsigned int> &vec);
    
protected:
    unsigned short int _rows;
    unsigned short int _cols;
    unsigned short int sparse_form_size;
    char* _graph_name;
    sparse_matrix_element<DT>** sparse_form;
    DT** _edges;
    
public:
    //Constructors
    DenseMatrix();
    DenseMatrix(const std::string &file_path);
    DenseMatrix(int rows, int cols);
    DenseMatrix(const DenseMatrix<DT>&);


    //Destructor
    virtual ~DenseMatrix();
    //Accessors
    const int getNumberOfRows();
    const int getNumberOfColumns();
    const bool isSquare();
    const bool isSymmetric();
    sparse_matrix_element<DT>** getSparseForm();
    int getSparseFormSize();
    DT getFrobNorm();
    const char* getGraphName();
    void setGraphName(const char* name);
    DenseMatrix<DT>* getScatteredSelection(std::vector<int>& vec_A, std::vector<int> vec_B);
    double* getTopEigenVector();
    vector<int>* getNeighbors(int vertex);
    //Mutators
    DenseMatrix<DT> transpose();
    DenseMatrix<DT>* kron(DenseMatrix<DT>& matrix);
    DT* sum_rows();
    DenseMatrix<DT>* vec_times_mat(DT*, int);
    DenseMatrix<DT>* mat_times_vec(DT*, int);
    
    //operators
    DT* operator[](int);
    friend std::ostream& operator<< <> (std::ostream& stream, const DenseMatrix<DT>& matrix);
    void operator= (const DenseMatrix<DT>&);
    DenseMatrix<DT> operator-(DenseMatrix<DT>& other_matrix);
    DenseMatrix<DT> operator*(DenseMatrix<DT>& other_matrix);
};

//==========================================================CONSTRUCTORS============================================================

template <typename DT>
DenseMatrix<DT>::DenseMatrix()
{
    this->_rows = DEFAULT_MATRIX_SIZE;
    this->_cols = DEFAULT_MATRIX_SIZE;
    if (!_initilalizeMatrix())
    {
        throw OutOfMemoryException();
    }
}

template<typename DT>
DenseMatrix<DT>::DenseMatrix(const std::string &file_path)
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
    file_reader >> tmp_x;      //skip the number of lines entry
    _initilalizeMatrix();
  
    
    while (!file_reader.eof())
    {
        file_reader >> tmp_x;
        file_reader >> tmp_y;
        this->_edges[tmp_x - 1][tmp_y - 1] = 1;
    }
    
    file_reader.close();
    
}

template <typename DT>
DenseMatrix<DT>::DenseMatrix(int rows, int cols)
{
    this->_rows = rows;
    this->_cols = cols;
    if (!_initilalizeMatrix())
    {
        throw OutOfMemoryException();
    }
}

template <typename DT>
DenseMatrix<DT>::DenseMatrix(const DenseMatrix<DT>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
template <typename DT>
DenseMatrix<DT>::~DenseMatrix()
{
    if (_edges != NULL)
    {
        for(int i=0; i < this->_rows; i++)
        {
            delete [] _edges[i];
        }
        delete [] _edges;
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
template <typename DT>
const int DenseMatrix<DT>::getNumberOfRows()
{
    return _rows;
}

template <typename DT>
const int DenseMatrix<DT>::getNumberOfColumns()
{
    return _cols;
}

template <typename DT>
const bool DenseMatrix<DT>::isSquare()
{
    return (this->_rows == this->_cols);
}

template <typename DT>
const bool DenseMatrix<DT>::isSymmetric()
{
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
template <typename T>
T DenseMatrix<T>::getFrobNorm(){
  T ret_val=0;

  for(int i=0;i<this->getNumberOfRows();i++){
    for(int j=0;j<this->getNumberOfColumns();j++){
      ret_val+= (*this)[i][j] * (*this)[i][j];
    }
  }

  return ret_val;
}

template <typename DT>
sparse_matrix_element<DT>** DenseMatrix<DT>::getSparseForm(){
    
    if(sparse_form==NULL){
        //   #pragma omp parallel for
        for(int i=0;i<_rows;i++){
            for(int j=0;j<_cols;j++){
                if(this->_edges[i][j]!=0){
                    //       #pragma omp critical
                    {
                        sparse_form_size++;
                    }
                }
            }
        }
        
        sparse_form= new sparse_matrix_element<DT>* [sparse_form_size];
        int counter=0;
        
        for(int i=0;i<_rows;i++){
            for(int j=0;j<_cols;j++){
                
                
                if(this->_edges[i][j]!=0){
                    
                    sparse_matrix_element<DT> *ele= new sparse_matrix_element<DT>;
                    ele->row_index=i;
                    ele->col_index=j;
                    ele->value=this->_edges[i][j];
                    sparse_form[counter]=ele;
                    counter++;
                    
                }
            }
        }
    }
    
    return sparse_form;
}

template <typename DT>
DenseMatrix<DT>* DenseMatrix<DT>::kron(DenseMatrix<DT>& matrix)
{
    // also checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
        // Throw something
    }
    
    //Initializing and allocating the product matrix
    int prod_size = this->_rows * matrix._rows;
    DenseMatrix<DT>* prod_matrix = new DenseMatrix<DT>(prod_size, prod_size);
    
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
                    (*prod_matrix)[(i_outer*matrix._rows) + i_inner][(j_outer*matrix._cols)+j_inner] =
                    this->_edges[i_outer][j_outer] * matrix._edges[i_inner][j_inner];
                }
            }
        }
    }
    
    return prod_matrix;
}

template <typename DT>
DenseMatrix<DT>* DenseMatrix<DT>::getScatteredSelection(std::vector<int>& vec_A, std::vector<int> vec_B)
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
    DenseMatrix<DT>* res_matrix = new DenseMatrix<DT>(num_in_A, num_in_B);

    int counter = 0;
    
    for (int i=0; i< vec_A.size(); i++)
    {
        for(int j=0; j< vec_B.size(); j++)
        {
            if ( vec_A[i] == 1 && vec_B[j] ==1)
            {
                (*res_matrix)[counter/num_in_B][counter%num_in_B] = this->_edges[i][j];
                counter++;
            }
        }
    }
    return res_matrix;
}


template <typename DT>
vector<int>* DenseMatrix<DT>::getNeighbors(int vertex){
  vector<int>* neighbors=new vector<int>();
  int counter=0;
  for(int i=0;i<this->getNumberOfRows();i++){
    if(this->_edges[i][vertex]==1){
      neighbors->push_back(i);
      //counter++;
    }
  }
 
  return neighbors;
} 

template <typename DT>
double* DenseMatrix<DT>::getTopEigenVector(){
#ifdef ARPACK
// 	std::cout<< "Using ARPACK" << std::endl;
    int arr_size= (this->_rows*(this->_rows+1))/2;
    double nzval[arr_size];
    int counter=0;
    
    for(int i=0; i<this->_rows; i++){
        for(int j=i;j < this->_rows ;j++){
            nzval[counter]= this->_edges[i][j];
            counter++;
        }
    }
    
    ARdsSymMatrix<double> ARMatrix(this->_rows,nzval,'L');
    ARluSymStdEig<double> eigProb(1, ARMatrix, "LM",10);
    eigProb.FindEigenvectors();
    
    double* eigenVec = new double[eigProb.GetN()];
     
    for (int i=0; i < eigProb.GetN() ; i++)
    {
        eigenVec[i] = eigProb.Eigenvector(0,i);
    }

    return eigenVec;
#endif
#ifdef EIGEN
// 		std::cout<< "Using EIGEN" << std::endl;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> A_eigen = Eigen::MatrixXd::Zero(this->_rows, this->_cols);
		for (long i=0; i< this->_rows; i++) 
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
		double* eigenVec = new double[this->_rows];

		for ( int i=0; i < this->_cols; i++)
		{	
			if (evals_eigen(i) == evals_eigen.maxCoeff())
			{
				for (int j= 0; j < this->_rows ; j++)
				{
					eigenVec[j] = evecs_eigen(j,i);
				}
				break;
			}
		}
	return eigenVec;
#endif
}

template <typename DT>
int DenseMatrix<DT>:: getSparseFormSize(){
    return sparse_form_size;
}

template <typename DT>
const char* DenseMatrix<DT>::getGraphName()
{
	return this->_graph_name;
}

//===========================================================MUTATORS================================================================

template <typename DT>
void DenseMatrix<DT>::setGraphName(const char* name)
{
	std::strcpy(this->_graph_name, name);
}
/*
 * returns an array of floats where each array element corresponds to sum of entries in a matrix row
 * @pram: pointer to 2-d array of floats
 * @pram: number of rows in graph
 * @pram: number of columns in graph
 */
template <typename DT>
DT* DenseMatrix<DT>::sum_rows(){
    DT* ret_arr= new DT[this->_rows];
    
#pragma omp parallel for
    for(int i=0;i<_rows;i++){
        ret_arr[i]=sum_array(this->_edges[i],_cols);
    }
    
    return ret_arr;
}


/*
 * returns the product of a vector and a matrix
 * @pram: an array of floats representing the vector
 * @pram: the size of the vector
 * @pram: a 2-d array of floats represeting the matrix
 * @pram: the number of rows in the matirx
 * @pram: the number of cols in the matrix
 */
template <typename DT>
DenseMatrix<DT>* DenseMatrix<DT>:: vec_times_mat(DT* vec, int vec_size){
    
    DenseMatrix<DT>* ret_matrix= new DenseMatrix<DT>(*this);
    if(_rows!=vec_size){
       fprintf(stderr,"dimension mismatch");
    }
    else{
        
        
#pragma omp for
        for(int i=0;i<_rows;i++){
			DT* tmp = ret_matrix->_edges[i];
            ret_matrix->_edges[i] = scalar_multiplication(tmp, _cols, vec[i]);
			delete [] tmp;
        }
    }
    return ret_matrix;
}

/*
 * returns the product of a matrix and a vector
 * @pram: pointer to the 2-d array representing the matrix
 * @pram: the number of rows in the matrix
 * @pram: the number of cols in the matrix
 * @pram: pointer to the array represting a vector
 * @pram: the size of the array
 */
template <typename DT>
DenseMatrix<DT>* DenseMatrix<DT>::mat_times_vec(DT* vec, int vec_size){
    
    DenseMatrix<DT>* ret_matrix=new DenseMatrix<DT>(*this);
    
    if(_cols!=vec_size){
        fprintf(stderr,"dimension mismatch");
    }
    else{
        
#pragma omp for
        for(int i=0;i<_cols;i++){
            for(int j=0;j<_rows;j++){
                ret_matrix->_edges[j][i]=this->_edges[j][i]*vec[i];
            }
        }
    }
    
    return ret_matrix;
}


template <typename DT>
DenseMatrix<DT> DenseMatrix<DT>::transpose(){
  DenseMatrix<DT> ret_matrix(this->getNumberOfColumns(),this->getNumberOfRows());


  for(int i=0;i<this->getNumberOfRows();i++){
    for(int j=0;j<this->getNumberOfColumns();j++){
      ret_matrix[j][i]=(*this)[i][j];
    }
  }

  return ret_matrix;
}


//==========================================================OPERATORS================================================================
template <typename DT>
std::ostream& operator<<(std::ostream& stream, const DenseMatrix<DT>& matrix)
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

template <typename DT>
DT* DenseMatrix<DT>::operator[](int index)
{
    return this->_edges[index];
}

template <typename DT>
DenseMatrix<DT> DenseMatrix<DT>::operator*(DenseMatrix<DT>& other_matrix){
  DenseMatrix<DT> ret_matrix(this->getNumberOfRows(),other_matrix.getNumberOfColumns());
  DT ret_val;
  for(int i=0;i<this->getNumberOfRows();i++){
    for(int j=0;j<other_matrix.getNumberOfColumns();j++){
      ret_val=0;
      for(int k=0;k<this->getNumberOfColumns();k++){
	ret_val+=(*this)[i][k]*other_matrix[k][j];
      }
      ret_matrix[i][j]=ret_val;
    }
  }
  
  return ret_matrix;
}



template <typename T>
DenseMatrix<T> DenseMatrix<T>::operator-(DenseMatrix<T>& other_matrix){
  DenseMatrix<T> ret_matrix(this->getNumberOfRows(),this->getNumberOfColumns());
  

  for(int i=0;i<this->getNumberOfRows();i++)
  {   
      for(int j=0;j<this->getNumberOfColumns();j++)
      {
	ret_matrix[i][j]=(*this)[i][j]-other_matrix[i][j];

      }
  }

  return ret_matrix;
}

template <typename DT>
void DenseMatrix<DT>::operator=(const DenseMatrix<DT>& matrix)
{
    _copy(matrix);
}

//===========================================================PRIVATE=================================================================
template <typename DT>
void DenseMatrix<DT>::_copy(const DenseMatrix<DT>& matrix)
{
    //delete this;
    this->_rows = matrix._rows;
    this->_cols = matrix._cols;
    if (!_initilalizeMatrix())
    {
        throw OutOfMemoryException();
    }
    for(int i=0; i< this->_rows; i++)
    {
        for(int j=0; j < this->_cols; j++)
        {
            this->_edges[i][j] = matrix._edges[i][j];
        }
    }
}

template <typename DT>
bool DenseMatrix<DT>::_initilalizeMatrix()
{
    this->sparse_form_size = 0;
    this->sparse_form = NULL;
    bool success = true;
    
    this->_edges = new DT*[this->_rows];
    success &= (this->_edges != NULL);               //checking for memory issues.
    
    for(int i=0; i < this->_rows; i++)
    {
        this->_edges[i] = new DT[this->_cols];
        for(int j=0; j< this->_cols; j++)
        {
            this->_edges[i][j] = 0;
        }
        
        success &= (this->_edges[i] != NULL);        //checking for memory issues.
    }
    
    return success;
}

template <typename DT>
void DenseMatrix<DT>::_split(const std::string &line, std::vector<unsigned int> &vec)
{
    vec.clear();
    std::istringstream iss(line);
    std::vector<string> tokens;
    std::copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));
    
    for (int i=0; i < tokens.size(); i++)
    {
        vec.push_back(std::atoi(tokens[i].c_str()));
    }
}

//===================================================================================================================================
       


#endif
    
