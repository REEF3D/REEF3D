/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"reinidisc.h"
#include"ddweno.h"
#include"vec.h"

class picard;
class cpt;

using namespace std;

#ifndef REINIDISC_F2_H_
#define REINIDISC_F2_H_

class reinidisc_f2 : public reinidisc, public ddweno
{
public:
	reinidisc_f2(lexer* p);
	virtual ~reinidisc_f2();
	virtual void start(lexer*, fdm*, ghostcell*, vec&, vec&,int);
	
private:
	void disc(lexer*, fdm*, ghostcell*, vec&, vec&, int*, cpt&);
	
	double xgrad,ygrad,zgrad;
    double xmin,xplus,ymin,yplus,zmin,zplus;
	double dxmin,dxplus,dymin,dyplus,dzmin,dzplus;
	double ux,uy,uz,ddt;
	double lsv,dv,S0;
	
	double dx, dy, dz, dnorm, lnorm, sign;
	double sx,sy,sz,snorm,op;
    double lsSig;
	
	double deltax,denom;
};

#endif
