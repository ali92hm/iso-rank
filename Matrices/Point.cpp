/*********************************************************************************
 * Implementation file for the Point class                                       *
 *********************************************************************************/
#include "Point.h"

//==========================================================CONSTRUCTORS============================================================
 /*
  * Default constructor: 
  * Initializes x and y to zero
  */
inline Point::Point()
{
    this->_x = 0;
    this->_y = 0;
}

 /*
  * Constructor with given x and y 
  * @pram int x
  * @pram int y
  */
inline Point::Point(int x,int y)
{
    this->_x = x;
    this->_y = y;
}

 /*
  * Copy constructor
  * @pram edge
  */
inline Point::Point(const Point& edge)
{
    _copy(edge);
}
//==========================================================DESTRUCTOR==============================================================

 /*
  * Destructor
  */
inline Point::~Point()
{

}
//===========================================================ACCESSORS===============================================================
 /*
  * X coordinate accessor
  * @retun int x
  */
inline int Point::getX() const
{
    return this->_x;
}

 /*
  * Y coordinate accessor
  * @retun int y
  */
inline int Point::getY() const
{
    return this->_y;
}

//===========================================================MUTATORS================================================================
 /*
  * X coordinate mutator
  * @pram int x
  */
inline void Point::setX(int x)
{
    this->_x = x;
}

 /*
  * Y coordinate mutator
  * @pram int y
  */
inline void Point::setY(int y)
{
    this->_y = y;
}

//==========================================================OPERATORS================================================================

/*
 * Overloaded > operator for point.
 * Points are sorted based their x value, y value is used if the points have the same x
 */
inline bool Point::operator>(const Point& rhs) const
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

/*
 * Overloaded < operator for point.
 * Points are sorted based their x value, y value is used if the points have the same x
 */
inline bool Point::operator<(const Point& rhs) const
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

/*
 * Overloaded equal operator for point.
 * Two pints are equal if they have the same x and y as each other
 */
inline bool Point::operator==(const Point& rhs) const
{
    return ((this->_x == rhs._x) && (this->_y == rhs._y));
}

/*
 * Overloaded assignment operator.
 */
inline void Point::operator=(const Point& rhs)
{
    _copy (rhs);
}

/*
 * overloaded print "<<" operator for point.
 */
inline std::ostream& operator<<(std::ostream& stream, const Point& p)
{
    stream<< "(" << p._x << ", " << p._y << ")" << std::endl;
    return stream;
}

//===========================================================PRIVATE=================================================================
 /*
  * Copy a point
  * @pram Point p
  */
inline void Point::_copy(const Point& p)
{
    this->_x = p._x;
    this->_y = p._y;
}

//===================================================================================================================================
