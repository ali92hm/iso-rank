/*********************************************************************************
 * This class represent an element in a sparse matrix by keeping track of its row *
 * and column number (i, j) and the value at that point                           *
 *********************************************************************************/

#ifndef _SparseElement_h
#define _SparseElement_h

#include <iostream>
#include "MatrixExceptions.h"

template <typename T>
class SparseElement;

template <typename T>
std::ostream& operator<<(std::ostream&, const SparseElement<T>&);

/*
 * SparseElement class definition and method declaration.
 */
template <typename T>
class SparseElement
{
private:
    void _copy(const SparseElement<T>&);
    
protected:
    int _i;
    int _j;
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
    int getI() const;
    int getJ() const;
    T getValue() const;
    
    /**********
     *MUTATORS*
     **********/
    void setI(int);
    void setJ(int);
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
/*
  * Default constructor: 
  * Initializes i, j and value to 0.
  */
template <typename T>
inline SparseElement<T>::SparseElement()
{
    this->_i = 0;
    this->_j = 0;
    this->_value = 0;
}

 /*
  * Constructor with given i and j. Initializes value to 0.
  * @pram int i
  * @pram int j
  */
template <typename T>
inline SparseElement<T>::SparseElement(int i,int j)
{
    this->_i = i;
    this->_j = j;
    this->_value = 0;
}

 /*
  * Constructor with given i, j and value.
  * @pram int i
  * @pram int j
  * @pram T value
  */
template <typename T>
inline SparseElement<T>::SparseElement(int i, int j, T value)
{
    this->_i = i;
    this->_j = j;
    this->_value = value;
}

 /*
  * Copy constructor
  * @pram edge
  */
template <typename T>
inline SparseElement<T>::SparseElement(const SparseElement<T>& edge)
{
    _copy(edge);
}
//==========================================================DESTRUCTOR==============================================================
 /*
  * Destructor
  */
template <typename T>
inline SparseElement<T>::~SparseElement()
{
    
}
//===========================================================ACCESSORS===============================================================
 /*
  * i coordinate accessor
  * @retun int i
  */
template <typename T>
inline int SparseElement<T>::getI() const
{
    return this->_i;
}

 /*
  * j coordinate accessor
  * @retun int j
  */
template <typename T>
inline int SparseElement<T>::getJ() const
{
    return this->_j;
}

 /*
  * Value accessor
  * @retun T value
  */
template <typename T>
inline T SparseElement<T>::getValue() const
{
    return this->_value;
}
//===========================================================MUTATORS================================================================
 /*
  * i coordinate mutator
  * @pram int i
  */
template <typename T>
inline void SparseElement<T>::setI(int i)
{
    this->_i = i;
}

 /*
  * j coordinate mutator
  * @pram int j
  */
template <typename T>
inline void SparseElement<T>::setJ(int j)
{
    this->_j = j;
}

 /*
  * Value mutator
  * @pram T value
  */
template <typename T>
inline void SparseElement<T>::setValue(T value)
{
    this->_value = value;
}
//==========================================================OPERATORS================================================================
/*
 * Overloaded > operator for SparseElement.
 * SparseElements are sorted based their i value, j value is used if the points have the same i
 */
template <typename T>
inline bool SparseElement<T>::operator>(const SparseElement<T>& rhs) const
{
    if (this->_i == rhs._i)
    {
        return (this->_j > rhs._j);
    }
    else
    {
        return (this->_i > rhs._i);
    }
}

/*
 * Overloaded < operator for SparseElement.
 * SparseElements are sorted based their i value, j value is used if the points have the same i
 */
template <typename T>
inline bool SparseElement<T>::operator<(const SparseElement<T>& rhs) const
{
    if (this->_i == rhs._i)
    {
        return (this->_j < rhs._j);
    }
    else
    {
        return (this->_i < rhs._i);
    }
}

/*
 * Overloaded equal operator for SparseElement.
 * Two SparseElements are equal if they have the same i and j as each other
 */
template <typename T>
inline bool SparseElement<T>::operator==(const SparseElement<T>& rhs) const
{
    return ((this->_i == rhs._i) && (this->_j == rhs._j));
}

/*
 * Overloaded assignment operator.
 */
template <typename T>
inline SparseElement<T>& SparseElement<T>::operator=(const SparseElement<T>& rhs)
{
    _copy (rhs);
}

/*
 * overloaded print "<<" operator for SparseElement.
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const SparseElement<T>& edge)
{
    stream<< edge._value << " at (" << edge._i << ", " << edge._j << ")" << std::endl;
    return stream;
}

//===========================================================PRIVATE=================================================================
 /*
  * Copy a SparseElement
  * @pram SparseElement<T> edge
  */
template <typename T>
inline void SparseElement<T>::_copy(const SparseElement<T>& edge)
{
    this->_i = edge._i;
    this->_j = edge._j;
    this->_value = edge._value;
}

//===================================================================================================================================
#endif
