/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"solver.h"
#include"increment.h"
#include"vec.h"

class cpt;

using namespace std;

#ifndef SIP_H_
#define SIP_H_

class sip : public solver, public increment
{
public:

	sip(lexer*,fdm*,ghostcell*);
	virtual ~sip();
	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, vec&, int, int, double);
    virtual void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int, int, double);
	virtual void solve(lexer*,fdm*, ghostcell*, vec&, vec&, int, int, int&, int, double, cpt&);
	virtual void setup(lexer*,fdm*, ghostcell*,int, cpt&);
	virtual void fillxvec1(lexer*,fdm*,field&);
    virtual void fillxvec2(lexer*,fdm*,field&);
    virtual void fillxvec3(lexer*,fdm*,field&);
    virtual void fillxvec4(lexer*,fdm*,field&);
	virtual void finalize(lexer*,fdm*,field&,vec&,int);
    virtual void gcpara_update(lexer*,vec&,ghostcell*);
	virtual void gcupdate(lexer*,fdm*,ghostcell*,vec&,int,int);

private:
	
	vec le,ls,lb,lp,un,uw,ut,res;
	
	const double epsi;
	int *sizeM;

	int numiter;
	double resi,y,residual,aii,alpha;
	double p1,p2,p3;
	int margin;
	int gcval_press;
	int qn,n,count,q;
	int ijk_n,ijk_s,ijk_e,ijk_w,ijk_t,ijk_b,ijk_p;

};

#endif

