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

#include"nhflow_pressure.h"
#include"increment.h"

class density;
class solver;
class nhflow_poisson;

using namespace std;

#ifndef NHFLOW_PJM_HS_H_
#define NHFLOW_PJM_HS_H_

class nhflow_pjm_hs : public nhflow_pressure, public increment
{

public:

	nhflow_pjm_hs(lexer* p, fdm_nhf*);
	virtual ~nhflow_pjm_hs();

	virtual void start(lexer*,fdm_nhf*,solver*,ghostcell*,ioflow*,double*,double*,double*,double);
	virtual void ucorr(lexer*p,fdm_nhf*,double*,double);
	virtual void vcorr(lexer*p,fdm_nhf*,double*,double);
	virtual void wcorr(lexer*p,fdm_nhf*,double*,double);
	virtual void upgrad(lexer*,fdm_nhf*,slice&,slice&);
	virtual void vpgrad(lexer*,fdm_nhf*,slice&,slice&);
    virtual void wpgrad(lexer*,fdm_nhf*,slice&,slice&);
    
	void rhs(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
	void vel_setup(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
    void bedbc(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);


private:
	double starttime,endtime;
	int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
	
	void debug(lexer*,fdm_nhf*);
    
    density *pd;
};



#endif

