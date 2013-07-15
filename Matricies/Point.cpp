//
//  Point.cpp
//  Graph Matching
//
//  Created by Ali Hajimirza on 7/10/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "Point.h"

//==========================================================CONSTRUCTORS============================================================
Point::Point()
{
    this->_x = -1;
    this->_y = -1;
}

Point::Point(int x,int y)
{
    this->_x = x;
    this->_y = y;
}

Point::Point(const Point& edge)
{
    _copy(edge);
}
//==========================================================DESTRUCTOR==============================================================
Point::~Point()
{
    
}
//===========================================================ACCESSORS===============================================================
int Point::getX() const
{
    return this->_x;
}

int Point::getY() const
{
    return this->_y;
}

//===========================================================MUTATORS================================================================
void Point::setX(int x)
{
    this->_x = x;
}

void Point::setY(int y)
{
    this->_y = y;
}

//==========================================================OPERATORS================================================================

bool Point::operator>(const Point& rhs) const
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

bool Point::operator<(const Point& rhs) const
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


bool Point::operator==(const Point& rhs) const
{
    return ((this->_x == rhs._x) && (this->_y == rhs._y));
}

void Point::operator=(const Point& rhs)
{
    _copy (rhs);
}

std::ostream& operator<<(std::ostream& stream, const Point& p)
{
    stream<< "(" << p._x << ", " << p._y << ")" << std::endl;
    return stream;
}

//===========================================================PRIVATE=================================================================

void Point::_copy(const Point& edge)
{
    this->_x = edge._x;
    this->_y = edge._y;
}

//===================================================================================================================================
