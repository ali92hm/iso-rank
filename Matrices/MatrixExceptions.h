//
//  MatrixExceptions.h
//  Graph Matching
//
//  Created by Ali Hajimirza on 7/5/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef __MatrixExceptions__
#define __MatrixExceptions__

#include <iostream>
#include <stdexcept>
#include <cstring>

class MatrixException : public std::exception{};
class OutOfMemoryException : public MatrixException {};
class IndexOutOfBoundsException : public MatrixException {};
class DimensionMismatchException : public MatrixException {};

class SymMatrixException : public MatrixException {};
class NotASymmetricMatrixException : public MatrixException {};
class NotASquareMatrixException : public MatrixException {};


class MatrixReaderException : public MatrixException {};
class FileDoesNotExistException : public std::exception
{
protected:
    char* _file_name;
public:
    explicit FileDoesNotExistException(const char* file_name)
    {
        std::strcpy(this->_file_name, file_name);
    }
    
    virtual const char* what() const throw()
    { 
    	return this->_file_name;
    }
};


#endif /* defined(__Graph_Matching__MatrixExceptions__) */
