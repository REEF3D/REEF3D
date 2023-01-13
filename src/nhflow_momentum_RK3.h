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

#include"nhflow_momentum.h"
#include"bcmom.h"


using namespace std;

#ifndef NHFLOW_MOMENTUM_RK3_H_
#define NHFLOW_MOMENTUM_RK3_H_

class nhflow_momentum_RK3 : public nhflow_momentum, public bcmom
{
public:
	nhflow_momentum_RK3(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_momentum_RK3();
    
	virtual void start(lexer*, fdm_nhf*, ghostcell*, ioflow*, convection*, diffusion*, nhflow_pressure*, solver*, nhflow*, nhflow_fsf*);


    double *UDIFF,*URK1,*URK2;
    double *VDIFF,*VRK1,*VRK2;
    double *WDIFF,*WRK1,*WRK2;

private:

	void irhs(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
	void jrhs(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
	void krhs(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
	
	int gcval_u, gcval_v, gcval_w;
	double starttime;
};

#endif
