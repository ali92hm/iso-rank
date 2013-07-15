//
//  SparseElement.h
//  Graph Matching
//
//  Created by Ali Hajimirza on 7/9/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _SparseElement_h
#define _SparseElement_h

#include <iostream>
#include "MatrixExceptions.h"

template <typename T>
class SparseElement;

template <typename T>
std::ostream& operator<<(std::ostream&, const SparseElement<T>&);

template <typename T>
class SparseElement
{
private:
    void _copy(const SparseElement<T>&);
    
protected:
    int _x;
    int _y;
    T _value;
    
public:
    /**************
     *Constructors*
     **************/
    SparseElement();
    SparseElement(int, int);
    SparseElement(int, int, T);
    SparseElement(const SparseElement<T>&);
    
    /************
     *Destructor*
     ************/
    virtual ~SparseElement();
    
    /***********
     *ACCESSORS*
     ***********/
    int getX() const;
    int getY() const;
    T getValue() const;
    
    /**********
     *MUTATORS*
     **********/
    void setX(int);
    void setY(int);
    void setValue(T);
    
    /**********
     *OPERATORS*
     **********/
    bool operator>(const SparseElement<T>&) const;
    bool operator<(const SparseElement<T>&) const;
    bool operator==(const SparseElement<T>&) const;
    SparseElement<T>& operator=(const SparseElement<T>&);
    friend std::ostream& operator<< <> (std::ostream&, const SparseElement<T>&);
};

//==========================================================CONSTRUCTORS============================================================
template <typename T>
SparseElement<T>::SparseElement()
{
    this->_x = -1;
    this->_y = -1;
    this->_value = NULL;
}

template <typename T>
SparseElement<T>::SparseElement(int x,int y)
{
    if ( x < 0 || y < 0 ) //subject to chage
        throw IndexOutOfBoundsException();
    this->_x = x;
    this->_y = y;
    this->_value = 0;
}

template <typename T>
SparseElement<T>::SparseElement(int x, int y, T value)
{
    if ( x < 0 || y < 0 ) //subject to chage
        throw IndexOutOfBoundsException();
    this->_x = x;
    this->_y = y;
    this->_value = value;
}

template <typename T>
SparseElement<T>::SparseElement(const SparseElement<T>& edge)
{
    _copy(edge);
}
//==========================================================DESTRUCTOR==============================================================
template <typename T>
SparseElement<T>::~SparseElement()
{
    
}
//===========================================================ACCESSORS===============================================================
template <typename T>
int SparseElement<T>::getX() const
{
    return this->_x;
}

template <typename T>
int SparseElement<T>::getY() const
{
    return this->_y;
}

template <typename T>
T SparseElement<T>::getValue() const
{
    return this->_value;
}
//===========================================================MUTATORS================================================================
template <typename T>
void SparseElement<T>::setX(int x)
{
    if ( x < 0 )//subject to chage
        throw IndexOutOfBoundsException();
    this->_x = x;
}

template <typename T>
void SparseElement<T>::setY(int y)
{
    if ( y < 0 )//subject to chage
        throw IndexOutOfBoundsException();
    this->_y = y;
}

template <typename T>
void SparseElement<T>::setValue(T value)
{
    this->_value = value;
}
//==========================================================OPERATORS================================================================
template <typename T>
bool SparseElement<T>::operator>(const SparseElement<T>& rhs) const
{
    if (this->_x == rhs._x)
    {
        return (this->_y > rhs._y);
    }
    else
    {
        return (this->_x > rhs._x);
    }
}

template <typename T>
bool SparseElement<T>::operator<(const SparseElement<T>& rhs) const
{
    if (this->_x == rhs._x)
    {
        return (this->_y < rhs._y);
    }
    else
    {
        return (this->_x < rhs._x);
    }
}

template <typename T>
bool SparseElement<T>::operator==(const SparseElement<T>& rhs) const
{
    return ((this->_x == rhs._x) && (this->_y == rhs._y));
}

template <typename T>
SparseElement<T>& SparseElement<T>::operator=(const SparseElement<T>& rhs)
{
    _copy (rhs);
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const SparseElement<T>& edge)
{
    stream<< edge._value << " at (" << edge._x << ", " << edge._y << ")" << std::endl;
    return stream;
}

//===========================================================PRIVATE=================================================================
template <typename T>
void SparseElement<T>::_copy(const SparseElement<T>& edge)
{
    this->_x = edge._x;
    this->_y = edge._y;
    this->_value = edge._value;
}

//===================================================================================================================================
#endif
