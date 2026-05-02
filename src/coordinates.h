/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef COORDINATES_H_
#define COORDINATES_H_

#include"increment.h"

class fdm;
class lexer;
class field;

using namespace std;

class coordinates : virtual public increment
{
public:
    coordinates(lexer*);
	virtual ~coordinates();
    
    // world to model
    double Xin(double,double);
    double Yin(double,double);
    
    // model to world
    double Xout(double,double);
    double Yout(double,double);
    
private:
    lexer *p;
    
    double pos;
    
    int stop,count;
    int ii,jj,kk;
    
    int is,ie,iloc;
    int js,je,jloc;
    int ks,ke,kloc;

};

#endif
