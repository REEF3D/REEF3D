/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"vec.h"
#include"solver.h"
#include"increment.h"


using namespace std;

#ifndef BICGSTAB_H_
#define BICGSTAB_H_


class bicgstab : public solver, public increment
{
public:
	bicgstab(lexer*,fdm*,ghostcell*,int);

	virtual ~bicgstab();

	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, vec&, int, int,double);
    virtual void startF(lexer*, fdm_fnpf*, ghostcell*, double*, vec&, matrix_diag&, int, int, double);
    
	virtual void solve(lexer*,fdm*, ghostcell*, vec&, vec&, int, int,int&,int,double, cpt&);
	virtual void setup(lexer*,fdm*, ghostcell*,int,cpt&);
	
	virtual void fillxvec1(lexer*,fdm*,field&);
    virtual void fillxvec2(lexer*,fdm*,field&);
    virtual void fillxvec3(lexer*,fdm*,field&);
    virtual void fillxvec4(lexer*,fdm*,field&);
	virtual void finalize(lexer*,fdm*,field&,vec&,int);

    virtual void gcpara_update(lexer*,vec&,ghostcell*);
	virtual void gcupdate(lexer*,fdm*,ghostcell*,vec&,int,int,int);

	virtual double res_calc(lexer*,fdm*, vec&, ghostcell*,cpt&);
	virtual void matvec_axb(lexer*,fdm*, vec&, vec&, cpt&);
	virtual void matvec_std(lexer*,fdm* a, vec&, vec&,cpt&);
	
	
	solver *precon,*precon123,*precon4;
	

private:

	vec sj,rj,r0,vj,tj,pj,precoeff,ph,sh;
	
	int *sizeM,*range;

	const double epsi;

	int count;
	int margin;
	
	double alpha,beta,w1,w2,w,residual,norm_vj,norm_r0,norm_sj,norm_rj ;
    double r_j1, r_j, sigma ;
};

#endif

