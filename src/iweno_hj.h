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

#ifndef IWENO_HJ_H_
#define IWENO_HJ_H_

#include"convection.h"
#include"increment.h"

class flux;

using namespace std;

class iweno_hj : public convection, public increment
{
public:
	iweno_hj (lexer*);
	virtual ~iweno_hj();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    void wenoloop1(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop2(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop3(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop4(lexer*,fdm*,field&,int,field&,field&,field&);

	void is_south(field&);
	void is_north(field&);
	void is_east(field&);
	void is_west(field&);
	void is_top(field&);
	void is_bottom(field&);

	void alpha_calc();

	void weight_calc();

	void aij(fdm*, field&,field&,field&,field&);
	void aij_south(lexer*,fdm*,field&, field&);
	void aij_north(lexer*,fdm*,field&, field&);
	void aij_east(lexer*,fdm*,field&, field&);
	void aij_west(lexer*,fdm*,field&, field&);
	void aij_top(lexer*,fdm*,field&, field&);
	void aij_bottom(lexer*,fdm*,field&, field&);

	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten,epsilon,deltin;


	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double umin, umax, uplus;
	int count,rocount,countN,coliN,aiipos;
	int *range;
    
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    double iadvec,jadvec,kadvec;
    
    flux *pflux;

};

#endif
