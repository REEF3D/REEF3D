/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"convection.h"
#include"weno_nug_func.h"

class flux;

#ifndef IWENO_HJ_NUG_H_
#define IWENO_HJ_NUG_H_

using namespace std;

class iweno_hj_nug : public convection, public weno_nug_func
{
public:
	iweno_hj_nug (lexer*);
	virtual ~iweno_hj_nug();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    void wenoloop1(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop2(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop3(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop4(lexer*,fdm*,field&,int,field&,field&,field&);

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
    
    double *DX,*DY,*DZ;

};

#endif
