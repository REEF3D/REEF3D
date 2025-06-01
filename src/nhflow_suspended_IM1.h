/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef NHFLOW_SUSPENDED_IM1_H_
#define NHFLOW_SUSPENDED_IM1_H_

#include"nhflow_suspended.h"
#include"increment.h"

using namespace std;

class nhflow_suspended_IM1 : public nhflow_suspended, public increment
{
public:
	nhflow_suspended_IM1(lexer *, fdm_nhf*);
	virtual ~nhflow_suspended_IM1();
	virtual void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, sediment_fdm*);
	virtual void ctimesave(lexer*, fdm_nhf*);
    
    void suspsource(lexer*,fdm_nhf*,double*,sediment_fdm*);
    void bcsusp_start(lexer*,fdm_nhf*,ghostcell*,sediment_fdm*,double*);
	void clearrhs(lexer*,fdm_nhf*);
    void fillconc(lexer*,fdm_nhf*,ghostcell*,sediment_fdm*);

	int gcval_susp;

	double *CONCN;

private:
    void timesource(lexer*, fdm_nhf*, double*);
    double starttime;
    void fill_wvel(lexer*,fdm_nhf*,ghostcell*,sediment_fdm*); 
    double *WVEL;
    
    int count,q;
};

#endif
