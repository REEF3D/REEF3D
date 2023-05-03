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

#include"solver2D.h"
#include"increment.h"

using namespace std;

#ifndef SFLOW_BICGSTAB_H_
#define SFLOW_BICGSTAB_H_

class sflow_bicgstab : public solver2D, public increment
{
public:

	sflow_bicgstab(lexer*,ghostcell*);
	virtual ~sflow_bicgstab();
	virtual void start(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int);
	virtual void solve(lexer*, ghostcell*, matrix2D&, vec2D&, vec2D&, int, int&);
    
    void fillxvec(lexer*,slice&,vec2D&);
	void finalize(lexer*,slice&);

	double res_calc(lexer*, matrix2D&, ghostcell*, double*);
	void matvec_axb(lexer*, matrix2D&, double*, double*);
	void matvec_std(lexer*, matrix2D&, double*, double*);
    
    void precon_setup(lexer*, matrix2D&,ghostcell*);
    void precon_solve(lexer*,ghostcell*,double*,double*);
	
    
private:

    double *sj,*rj,*r0,*vj,*tj,*pj,*ph,*sh,*x,*rhs,*aii;

    int num_iterations;
    double final_res_norm;
	int stencil_indices[7];
	int nentries;
   
	int numiter,count,q;
    
    
    //cg
	int *sizeS,*range;

	const double epsi;

	int margin;
    int ulast,vlast,wlast;
    int *flagslice;
	
	double alpha,beta,w1,w2,w,residual,norm_vj,norm_r0,norm_sj,norm_rj ;
    double r_j1, r_j, sigma;


};

#endif


