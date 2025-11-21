/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
/*
#ifndef SFLOW_PJM_HS_H_
#define SFLOW_PJM_HS_H_

#include"sflow_pressure.h"
#include"sflow_gradient.h"

class density;
class solver;
class patchBC_interface;

using namespace std;

class sflow_pjm_hs : public sflow_pressure, public sflow_gradient
{

public:

	sflow_pjm_hs(lexer* p, fdm2D*,patchBC_interface*);
	virtual ~sflow_pjm_hs();

	virtual void start(lexer*,fdm2D*,solver*,ghostcell*,ioflow*,slice&,double*,double*,double*,double);
	virtual void ucorr(lexer*p,fdm2D*,slice&,double*,double*,double);
	virtual void vcorr(lexer*p,fdm2D*,slice&,double*,double*,double);
	virtual void wcorr(lexer*p,fdm2D*,slice&,double*,double*,double);
	virtual void upgrad(lexer*,fdm2D*,slice&);
	virtual void vpgrad(lexer*,fdm2D*,slice&);
    virtual void wpgrad(lexer*,fdm2D*,slice&);
    
	void rhs(lexer*,fdm2D*,ghostcell*,double*,double*,double*,double);
	void vel_setup(lexer*,fdm2D*,ghostcell*,double*,double*,double*,double);
    void bedbc(lexer*,fdm2D*,ghostcell*,double*,double*,double*,double);


private:

	double starttime,endtime;
	int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
    double val, denom;
    
    
    density *pd;
    patchBC_interface *pBC;
};



#endif

*/