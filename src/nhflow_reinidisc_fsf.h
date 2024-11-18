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

#include"ddweno_nug_sig.h"

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

#ifndef NHFLOW_REINIDISC_FSF_H_
#define NHFLOW_REINIDISC_FSF_H_

class nhflow_reinidisc_fsf : public increment, public ddweno_nug_sig
{
public:
	nhflow_reinidisc_fsf(lexer* p);
	virtual ~nhflow_reinidisc_fsf();
    
	virtual void start(lexer*, fdm_nhf*, ghostcell*, double*, double*);
	
private:

	void disc(lexer*, fdm_nhf*, ghostcell*, double*, double*);
	
	double xmin,xplus,ymin,yplus,zmin,zplus;
	double dxmin,dxplus,dymin,dyplus,dzmin,dzplus;
	double uwx,uwy,uwz,ddt;
	double lsv,dv,lsSig;
	
	double dx, dy, dz, dnorm, sign;
	double sx,sy,sz,snorm,op;
	
	double deltax,denom;
};

#endif
