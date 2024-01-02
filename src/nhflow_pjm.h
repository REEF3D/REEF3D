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

#include"nhflow_pressure.h"
#include"increment.h"

class density;
class solver;
class nhflow_poisson;
class patchBC_interface;

using namespace std;

#ifndef NHFLOW_PJM_H_
#define NHFLOW_PJM_H_

class nhflow_pjm : public nhflow_pressure, public increment
{

public:

	nhflow_pjm(lexer*, fdm_nhf*, ghostcell*,patchBC_interface*);
	virtual ~nhflow_pjm();

	virtual void start(lexer*,fdm_nhf*,solver*,ghostcell*,ioflow*,slice&,double*,double*,double*,double);
	virtual void ucorr(lexer*p,fdm_nhf*,slice&,double*,double*,double);
	virtual void vcorr(lexer*p,fdm_nhf*,slice&,double*,double*,double);
	virtual void wcorr(lexer*p,fdm_nhf*,slice&,double*,double*,double);
	virtual void upgrad(lexer*,fdm_nhf*,slice&);
	virtual void vpgrad(lexer*,fdm_nhf*,slice&);
    virtual void wpgrad(lexer*,fdm_nhf*,slice&);
    
    void rhs(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
	void vel_setup(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);
    void bedbc(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double);

private:
    
	double starttime,endtime;
    const double teta;
    int check;
	int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
    int solver_id;
    double val, denom;

    density *pd;
    nhflow_poisson *ppois;
    patchBC_interface *pBC;

};



#endif

