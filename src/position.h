/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"increment.h"

class fdm;
class lexer;
class field;

#ifndef POSITION_H_
#define POSITION_H_

using namespace std;

class position : virtual public increment
{
public:
    position(lexer*);
	virtual ~position();
    
    // xyz
    double pos_x();
	double pos_y();
	double pos_z();
	
	double pos1_x();
	double pos1_y();
	double pos1_z();
	
	double pos2_x();
	double pos2_y();
	double pos2_z();
	
	double pos3_x();
	double pos3_y();
	double pos3_z();
	
	double posnode_x();
	double posnode_y();
	double posnode_z();
    
    // ijk
    int posf_i(double);
    int posf_j(double);
    int posf_k(double);
    
    int posc_i(double);
    int posc_j(double);
    int posc_k(double);
    
    int ihalf(int,int);
    
    int conv(double);
    
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
