/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"convection.h"
#include"increment.h"

class flux;

#ifndef HRIC_MOD_H_
#define HRIC_MOD_H_

using namespace std;

class hric_mod : public convection,  public increment
{

public:

	hric_mod (lexer *);
	virtual ~hric_mod();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&);
	
	double cface(lexer*,fdm*,field&,int,int,double);

	double dx,dy,dz;
    double fx2,fy2,fz2;
    double fx1,fy1,fz1;
	double ul,ur,vl,vr,wl,wr;
	double L,phi;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    flux *pflux;

};

#endif
