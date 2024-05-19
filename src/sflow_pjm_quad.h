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

#include"sflow_pressure.h"
#include"increment.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class sflow_weno_hj;
class sflow_gradient_weno;

using namespace std;

#ifndef SFLOW_PJM_QUAD_H_
#define SFLOW_PJM_QUAD_H_

class sflow_pjm_quad : public sflow_pressure, public increment
{
public:
    sflow_pjm_quad(lexer*, fdm2D*,patchBC_interface*);
	virtual ~sflow_pjm_quad();
    
	virtual void start(lexer*, fdm2D*, ghostcell*, solver2D*, ioflow*, slice&, slice&, slice&, slice&, slice&, slice&, double);
	virtual void upgrad(lexer*, fdm2D*, slice&, slice&);
	virtual void vpgrad(lexer*, fdm2D*, slice&, slice&);
    virtual void wpgrad(lexer*, fdm2D*, slice&, slice&);
    
    virtual void ucorr(lexer*,fdm2D*,slice&,slice&,double);
	virtual void vcorr(lexer*,fdm2D*,slice&,slice&,double);
	virtual void wcorr(lexer*,fdm2D*,double,slice&,slice&,slice&);
    virtual void wcalc(lexer*,fdm2D*,double,slice&,slice&,slice&);
    
    void rhs(lexer*, fdm2D*, slice&, slice&, slice&, double);
    
    void poisson(lexer*,fdm2D*,double);
    
private:
    void quad_prep(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,double);
    void quad_calc(lexer*,fdm2D*,slice&,slice&,slice&,slice&,double);
    
    
	double starttime,endtime;
    int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
    
    double sqd;
	double theta;
	double solvtime,ptime;
    double wd_criterion;
    
    slice4 phi4,press_n;
	sflow_weno_hj *disc;
    sflow_gradient_weno *pgrad;
    patchBC_interface *pBC;
    
    slice1 Ps;
    slice2 Qs;

};

#endif
