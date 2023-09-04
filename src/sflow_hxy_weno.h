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

#include"sflow_hxy_disc.h"
#include"increment.h"

class sflow_flux;

#ifndef SFLOW_HXY_WENO_H_
#define SFLOW_HXY_WENO_H_

using namespace std;

class sflow_hxy_weno : public sflow_hxy_disc, public increment
{
public:
	sflow_hxy_weno(lexer*,patchBC_interface*);
	virtual ~sflow_hxy_weno();

	virtual void start(lexer*,slice&,slice&,slice&,int*,slice&,slice&,slice&);

private:

	double fx(lexer*, slice&, int, double);
	double fy(lexer*, slice&, int, double);
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
    int qq;
    
    patchBC_interface *pBC;

};

#endif
