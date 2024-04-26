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

using namespace std;

#ifndef SFLOW_IWENO_HJ_H_
#define SFLOW_IWENO_HJ_H_

class sflow_iweno_hj : public sflow_convection, public increment
{
public:
	sflow_iweno_hj(lexer*);
	virtual ~sflow_iweno_hj();

	virtual void start(lexer*,fdm2D*,slice&,int,slice&,slice&);

private:
    double aij(lexer*, fdm2D*, slice&, int, slice&, slice&);
    
    void wenoloop1(lexer*,fdm2D*,slice&,int,slice&,slice&);
    void wenoloop2(lexer*,fdm2D*,slice&,int,slice&,slice&);
    void wenoloop4(lexer*,fdm2D*,slice&,int,slice&,slice&);

	void iqmin(lexer*, fdm2D*, slice&, int);
	void jqmin(lexer*, fdm2D*, slice&, int);
	void iqmax(lexer*, fdm2D*, slice&, int);
	void jqmax(lexer*, fdm2D*, slice&, int);
    
	void aij_south(lexer*,fdm2D*,slice&, slice&);
	void aij_north(lexer*,fdm2D*,slice&, slice&);
	void aij_east(lexer*,fdm2D*,slice&, slice&);
	void aij_west(lexer*,fdm2D*,slice&, slice&);
    
	void is_south(slice&);
	void is_north(slice&);
	void is_east(slice&);
	void is_west(slice&);

    void alpha_calc();

	void weight_calc();


	double L,grad;
	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon;
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
    
    const double deltin;
    int count;
};

#endif
