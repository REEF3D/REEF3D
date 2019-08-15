/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"vec.h"
#include"solver.h"
#include"increment.h"
#include"cpt.h"

class grid;

using namespace std;

#ifndef BICGSTAB_WIDE_H_
#define BICGSTAB_WIDE_H_


class bicgstab_wide : public solver, public increment
{
public:
	bicgstab_wide(lexer*);

	virtual ~bicgstab_wide();

	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, vec&, int, int,double);
    virtual void startF(lexer*, fdm_fnpf*, ghostcell*, double*, vec&, matrix_diag&, int, int, double);
    virtual void solveF(lexer*, fdm_fnpf*, ghostcell*, vec&, vec&, matrix_diag&, int, int, double);
    
	virtual void solve(lexer*,fdm*, ghostcell*, vec&, vec&, int, int,int&,int,double, cpt&);
	virtual void setup(lexer*,fdm*, ghostcell*,int,cpt&);


    
private:
    
    void fillxvec(lexer*,double*,vec&);
	void finalize(lexer*,double*,vec&);

    void gcpara_update(lexer*,vec&,ghostcell*);
	void gcupdate(lexer*,ghostcell*,vec&,int,int,int);

	double res_calc(lexer*, matrix_diag&, vec&, vec&, ghostcell*,cpt&);
	void matvec_axb(lexer*, matrix_diag&, vec&, vec&, vec&, cpt&);
	void matvec_std(lexer*, matrix_diag&, vec&, vec&, cpt&);
    
    void precon(lexer*, vec&, vec&,cpt&);
    void precon_setup(lexer*, matrix_diag&,cpt&);
	

private:

	vec xvec,sj,rj,r0,vj,tj,pj,precoeff,ph,sh,aii;
	
	int *sizeM,*range;

	const double epsi;

	int count;
	int margin;
    int solveriter;
	
	double alpha,beta,w1,w2,w,residual,norm_vj,norm_r0,norm_sj,norm_rj ;
    double r_j1, r_j, sigma;
    
    grid *pgrid;
    cpt C4;
};

#endif

