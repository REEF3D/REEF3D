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

#include"reini.h"
#include"gradient.h"
#include"fieldint4.h"

class picard;

using namespace std;

#ifndef REINIVC_RK3_H_
#define REINIVC_RK3_H_

class reinivc_RK3 : public reini, gradient
{
public:
	reinivc_RK3(lexer* p);
	virtual ~reinivc_RK3();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);
    virtual void startV(fdm*,lexer*,vec&,ghostcell*,ioflow*);

	double dx, dy, dz, dnorm, sign;
	double sx,sy,sz,snorm,op;

private:
    picard *ppicard;
	fieldint4 iflag;

	void step(fdm*, lexer*);
	double dt;
	void disc(lexer *p, fdm*, field&);
	
	void interface_cells(lexer *p, fdm*, field&);
	void correction(lexer *p, fdm*, field&);
	void finalize(lexer *p, fdm*);

	double xmin,xplus,ymin,yplus,zmin,zplus;
	double dxmin,dxplus,dymin,dyplus,dzmin,dzplus;
	double uwx,uwy,uwz,ddt;
	double lsv,dv,lsSig;
	double starttime,endtime;

	int gcval_phi,gcval_ro,gcval_iniphi,reiniter;
	const double epsi;
    double deltax;
	
	int count,numvert,n;
	double dV1,dV,H,H0,eta;
	
	double *ls,*ls1,*ls0;
	int **ijk;
};

#endif
