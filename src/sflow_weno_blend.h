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

#ifndef SFLOW_WENO_BLEND_H_
#define SFLOW_WENO_BLEND_H_

#include"sflow_convection.h"
#include"increment.h"

class sflow_flux;

using namespace std;

class sflow_weno_blend : public sflow_convection, public increment
{
public:
	sflow_weno_blend(lexer*);
	virtual ~sflow_weno_blend();

	virtual void start(lexer*,fdm2D*,slice&,int,slice&,slice&);

private:
    double aij_flux(lexer*, fdm2D*, slice&, int, slice&, slice&);
    double aij_hj(lexer*, fdm2D*, slice&, int, slice&, slice&);

	virtual double fx_flux(lexer*, fdm2D*, slice&, int, double);
	virtual double fy_flux(lexer*, fdm2D*, slice&, int, double);
    virtual double fx_hj(lexer*, fdm2D*, slice&, int, double);
	virtual double fy_hj(lexer*, fdm2D*, slice&, int, double);
    
	void iqmin_flux(lexer*, fdm2D*, slice&, int);
	void jqmin_flux(lexer*, fdm2D*, slice&, int);
	void iqmax_flux(lexer*, fdm2D*, slice&, int);
	void jqmax_flux(lexer*, fdm2D*, slice&, int);
    
    void iqmin_hj(lexer*, fdm2D*, slice&, int);
	void jqmin_hj(lexer*, fdm2D*, slice&, int);
	void iqmax_hj(lexer*, fdm2D*, slice&, int);
	void jqmax_hj(lexer*, fdm2D*, slice&, int);


	double L,grad;
	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum;
	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double q1,q2,q3,q4,q5;
	double gradx, grady, gradz;
	double fu1,fv1,fu2,fv2;
    double iadvec,jadvec;

	void is(slice&);
	void alpha();
	void weight();
    
    sflow_flux *pflux,*pflux_hj;
    double ivel1,ivel2,jvel1,jvel2;
};

#endif
