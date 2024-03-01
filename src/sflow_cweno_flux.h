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

#include"sflow_convection.h"
#include"increment.h"

class sflow_flux;
class fdm2D;

using namespace std;

#ifndef SFLOW_CWENO_FLUX_H_
#define SFLOW_CWENO_FLUX_H_

class sflow_cweno_flux : public sflow_convection, public increment
{
public:
	sflow_cweno_flux(lexer*,fdm2D*);
	virtual ~sflow_cweno_flux();

	virtual void start(lexer*,fdm2D*,slice&,int,slice&,slice&);

private:
    double aij(lexer*, fdm2D*, slice&, int, slice&, slice&);

	virtual double fx(lexer*, fdm2D*, slice&, int, double);
	virtual double fy(lexer*, fdm2D*, slice&, int, double);
	void iqmin(lexer*, fdm2D*, slice&, int);
	void jqmin(lexer*, fdm2D*, slice&, int);
	void iqmax(lexer*, fdm2D*, slice&, int);
	void jqmax(lexer*, fdm2D*, slice&, int);


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
    double dx,dy;


	void is_min(fdm2D*,slice&);
    void is_max(fdm2D*,slice&);
	void alpha();
	void weight();
    
    sflow_flux *pflux;
    double ivel1,ivel2,jvel1,jvel2;
};

#endif
