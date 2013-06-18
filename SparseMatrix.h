//
//  SparseMatrix.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _SparseMatrix_h
#define _SparseMatrix_h

#define DEFAULT_MATRIX_SIZE 20

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
#include "util.h"
#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"


template <typename DT>
class SparseMatrix;

template <typename DT>
struct sparse_matrix_element{
    int row_index;
    int col_index;
    DT value;
};

class SparseMatrixException : std::exception{};
class OutOfMemoryException : SparseMatrixException {};
class NotASymmetricMatrixException : SparseMatrixException {};
class NotASquareMatrixException : SparseMatrixException {};
class MatrixReaderException : SparseMatrixException {};
class FileDoesNotExistException : MatrixReaderException {};
class CurruptedFileException : MatrixReaderException {};
class WrongFormatException : MatrixReaderException {};

template <typename DT>
std::ostream& operator<< (std::ostream&, const SparseMatrix<DT>&);

template <typename DT>
class SparseMatrix
{
private:
    std::ifstream _file_reader;
    void _copy(const SparseMatrix<DT>&);
    bool _initilalizeMatrix();
    static void _split(const std::string &line, std::vector<unsigned int> &vec);
    
protected:
    unsigned short int _rows;
    unsigned short int _cols;
    unsigned short int sparse_form_size;
    sparse_matrix_element<DT>** sparse_form;
    DT** _edges;
    
public:
    //Constructors
    SparseMatrix();
    SparseMatrix(std::string &file_path);
    SparseMatrix(int rows, int cols);
    SparseMatrix(const SparseMatrix<DT>&);
    //Destructor
    virtual ~SparseMatrix();
    //Accessors
    const int getNumberOfRows();
    const int getNumberOfColumns();
    const bool isSquare();
    const bool isSymmetric();
    sparse_matrix_element<DT>** getSparseForm();
    int getSparseFormSize();
    SparseMatrix<DT>* getScatteredSelection(std::vector<int>& vec_A, std::vector<int> vec_B);
    ARdsSymMatrix<DT>* SparseMatrix<DT>::getARMatrix();
    //Mutators
    SparseMatrix<DT>* kron(SparseMatrix<DT>& matrix);
    DT* sum_rows();
    SparseMatrix<DT>* vec_times_mat(DT*, int);
    SparseMatrix<DT>* mat_times_vec(DT*, int);
    
    //operators
    DT* operator[](int);
    friend std::ostream& operator<< <> (std::ostream& stream, const SparseMatrix<DT>& matrix);
    void operator= (const SparseMatrix<DT>&);
};

//==========================================================CONSTRUCTORS============================================================

template <typename DT>
SparseMatrix<DT>::SparseMatrix()
{
    this->_rows = DEFAULT_MATRIX_SIZE;
    this->_cols = DEFAULT_MATRIX_SIZE;
    if (!_initilalizeMatrix())
    {
        throw OutOfMemoryException();
    }
}

template<typename DT>
SparseMatrix<DT>::SparseMatrix(std::string &file_path)
{
    _file_reader.open(file_path.c_str());
    
    //Throw exception if file doesn't exists.
    if(_file_reader.bad())
    {
        _file_reader.close();
        throw FileDoesNotExistException();
    }
    
    //Throw exception if the input file has format errors
    if (_file_reader.fail())
    {
        _file_reader.close();
        throw CurruptedFileException();
    }
    
    std::string line;
    int size;
    std::vector<unsigned int> values(3);
    
    if (_file_reader.is_open())
    {
        if (_file_reader.good())
        {
            getline(_file_reader, line);
            _split(line, values);
            
            this->_rows = values[0];
            this->_cols = values[1];
            size = values[2];
        
            
            if (!_initilalizeMatrix())
            {
                _file_reader.close();
                throw OutOfMemoryException();
            }
        }
        else
        {
            delete this;
            _file_reader.close();
            throw CurruptedFileException();
        }
        int i = 0;
        while( _file_reader.good())
        {
            getline(_file_reader, line );
            if (line == "")
            {
                continue;
            }
            _split(line, values);
            this->_edges[values[0]-1][values[1]-1] = 1;
            i++;
        }
        if (size != i)
        {
            delete this;
            _file_reader.close();
            throw CurruptedFileException();
        }
    }
    else
    {
        delete this;
        _file_reader.close();
        throw MatrixReaderException();
    }
    
    _file_reader.close();
    
}

template <typename DT>
SparseMatrix<DT>::SparseMatrix(int rows, int cols)
{
    this->_rows = rows;
    this->_cols = cols;
    if (!_initilalizeMatrix())
    {
        throw OutOfMemoryException();
    }
}

template <typename DT>
SparseMatrix<DT>::SparseMatrix(const SparseMatrix<DT>& matrix)
{
    _copy(matrix);
}

//==========================================================DESTRUCTOR==============================================================
template <typename DT>
SparseMatrix<DT>::~SparseMatrix()
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
const int SparseMatrix<DT>::getNumberOfRows()
{
    return _rows;
}

template <typename DT>
const int SparseMatrix<DT>::getNumberOfColumns()
{
    return _cols;
}

template <typename DT>
const bool SparseMatrix<DT>::isSquare()
{
    return (this->_rows == this->_cols);
}

template <typename DT>
const bool SparseMatrix<DT>::isSymmetric()
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

template <typename DT>
sparse_matrix_element<DT>** SparseMatrix<DT>::getSparseForm(){
    
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
SparseMatrix<DT>* SparseMatrix<DT>::kron(SparseMatrix<DT>& matrix)
{
    // also checking for matrices to be square
    if (!this->isSquare() || !matrix.isSquare())
    {
        // Throw something
    }
    
    //Initializing and allocating the product matrix
    int prod_size = this->_rows * matrix._rows;
    SparseMatrix<DT>* prod_matrix = new SparseMatrix<DT>(prod_size, prod_size);
    
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
SparseMatrix<DT>* SparseMatrix<DT>::getScatteredSelection(std::vector<int>& vec_A, std::vector<int> vec_B)
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
    SparseMatrix<DT>* res_matrix = new SparseMatrix<DT>(num_in_A, num_in_B);

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
ARdsSymMatrix<DT>* SparseMatrix<DT>::getARMatrix()
{
    int sparse_size=this->_rows;
    
    int arr_size= sparse_size*(sparse_size+1)/2;
    DT nzval[arr_size];
    int counter=0;
    
    for(int i=0;i<sparse_size;i++){
        for(int j=0;j<=i;j++){
            nzval[counter]= this->_edges[i][j];
            counter++;
        }
    }
    
    ARdsSymMatrix<DT> *matrix = new ARdsSymMatrix<DT>(sparse_size, nzval,'L');
    return matrix;
}

template <typename DT>
int SparseMatrix<DT>:: getSparseFormSize(){
    return sparse_form_size;
}

//===========================================================MUTATORS================================================================

/*
 * returns an array of floats where each array element corresponds to sum of entries in a matrix row
 * @pram: pointer to 2-d array of floats
 * @pram: number of rows in graph
 * @pram: number of columns in graph
 */
template <typename DT>
DT* SparseMatrix<DT>::sum_rows(){
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
SparseMatrix<DT>* SparseMatrix<DT>:: vec_times_mat(DT* vec, int vec_size){
    
    SparseMatrix<DT>* ret_matrix= new SparseMatrix<DT>(*this);
    if(_rows!=vec_size){
       fprintf(stderr,"dimension mismatch");
    }
    else{
        
        
#pragma omp for
        for(int i=0;i<_rows;i++){
            ret_matrix->_edges[i] = scalar_multiplication(ret_matrix->_edges[i], _cols, vec[i]);
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
SparseMatrix<DT>* SparseMatrix<DT>::mat_times_vec(DT* vec, int vec_size){
    
    SparseMatrix<DT>* ret_matrix=new SparseMatrix<DT>(*this);
    
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

//==========================================================OPERATORS================================================================
template <typename DT>
std::ostream& operator<<(std::ostream& stream, const SparseMatrix<DT>& matrix)
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
DT* SparseMatrix<DT>::operator[](int index)
{
    return this->_edges[index];
}

template <typename DT>
void SparseMatrix<DT>::operator=(const SparseMatrix<DT>& matrix)
{
    _copy(matrix);
}

//===========================================================PRIVATE=================================================================
template <typename DT>
void SparseMatrix<DT>::_copy(const SparseMatrix<DT>& matrix)
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
bool SparseMatrix<DT>::_initilalizeMatrix()
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
void SparseMatrix<DT>::_split(const std::string &line, std::vector<unsigned int> &vec)
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
    
