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
#include"ddweno_nug.h"
#include"vec.h"

class picard;
class cpt;

using namespace std;

#ifndef REINIDISC_FSF_RIG_H_
#define REINIDISC_FSF_RIG_H_

class reinidisc_fsf_rig : public reinidisc, public ddweno_nug
{
public:
	reinidisc_fsf_rig(lexer* p);
	virtual ~reinidisc_fsf_rig();
	virtual void start(lexer*, fdm*, ghostcell*, vec&, vec&,int);
	
private:

	void disc(lexer*, fdm*, ghostcell*, vec&, vec&, int*, int, cpt&);
	
	double xmin,xplus,ymin,yplus,zmin,zplus;
	double dxmin,dxplus,dymin,dyplus,dzmin,dzplus;
	double uwx,uwy,uwz,ddt;
	double lsv,dv,lsSig;
	
	double dx, dy, dz, dnorm, sign;
	double sx,sy,sz,snorm,op;
	
	double deltax,denom;
};

#endif
