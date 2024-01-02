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

#include"density.h"
#include"increment.h"

class fdm;
class lexer;

#ifndef DENSITY_DF_H_
#define DENSITY_DF_H_


using namespace std;

class density_df : public density, virtual public increment
{

public:
    density_df(lexer*);
	virtual ~density_df();

	virtual double roface(lexer*,fdm*,int,int,int);
	
	double H,H_fb,roval,phival,fbval;
	int ii,jj,kk;
	const double epsi,eps;
    double chi;
    double r,s;

};

#endif




