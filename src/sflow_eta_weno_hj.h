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

#ifndef SFLOW_ETA_WENO_HJ_H_
#define SFLOW_ETA_WENO_HJ_H_

#include"sflow_eta_disc.h"
#include"increment.h"

class sflow_flux;

using namespace std;

class sflow_eta_weno_hj : public sflow_eta_disc, public increment
{
public:
	sflow_eta_weno_hj(lexer*);
	virtual ~sflow_eta_weno_hj();

	virtual void start(lexer*,slice&,int,slice&,slice&,slice&,slice&);

private:
    double aij(lexer*, slice&, int, slice&, slice& ,slice&);

	virtual double fx(lexer*, slice&, int, double);
	virtual double fy(lexer*, slice&, int, double);
	void iqmin(lexer*, slice&, int);
	void jqmin(lexer*, slice&, int);
	void iqmax(lexer*, slice&, int);
	void jqmax(lexer*, slice&, int);


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

	void is(slice&);
	void alpha();
	void weight();
    
    sflow_flux *pflux;
    double ivel1,ivel2,jvel1,jvel2;
    double iadvec,jadvec;
};

#endif
