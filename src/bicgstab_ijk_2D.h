/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"solver.h"
#include"increment.h"

using namespace std;

#ifndef BICGSTAB_IJK_2D_H_
#define BICGSTAB_IJK_2D_H_


class bicgstab_ijk_2D : public solver, public increment
{
public:
	bicgstab_ijk_2D(lexer*,fdm*,ghostcell*);

	virtual ~bicgstab_ijk_2D();

	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, int);
    virtual void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int);
    virtual void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    virtual void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    virtual void startM(lexer*, ghostcell*, double*, double*, double*, int);
    
	virtual void solve(lexer*, ghostcell*, vec&, matrix_diag&, int, int&,int,double);
	
	void fillxvec(lexer*,fdm*,field&,vec&);
	void finalize(lexer*,fdm*,field&);

	void fillxvecV(lexer*,double*,vec&);
	void finalizeV(lexer*,double*);

	double res_calc(lexer*,ghostcell*, double*, matrix_diag&);
	void matvec_axb(lexer*, double*, double*, matrix_diag&);
	void matvec_std(lexer*, double*, double*, matrix_diag&);
    
    void precon_setup(lexer*,ghostcell*, matrix_diag&);
    void precon_solve(lexer*,ghostcell*,double*,double*, matrix_diag&);
	

private:

	double *sj,*rj,*r0,*vj,*tj,*pj,*ph,*sh,*x,*rhs,*aii;
	
	int *sizeM,*range;

	const double epsi;

	int count;
	int margin;
    int ulast,vlast,wlast;
    int *flag;
    double stop_crit;
	
	double alpha,beta,w1,w2,w,residual,norm_vj,norm_r0,norm_sj,norm_rj ;
    double r_j1, r_j, sigma ;
    
};

#endif

