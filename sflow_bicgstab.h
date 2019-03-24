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

#define HYPRE_COMPILATION
#ifdef  HYPRE_COMPILATION

#include"solver2D.h"
#include"increment.h"
#include"vec2D.h"

using namespace std;

#ifndef SFLOW_BICGSTAB_H_
#define SFLOW_BICGSTAB_H_

class sflow_bicgstab : public solver2D, public increment
{
public:

	sflow_bicgstab(lexer*,fdm2D*,ghostcell*);
	virtual ~sflow_bicgstab();
	virtual void start(lexer*,fdm2D*, ghostcell*, slice&, vec2D&, vec2D&, int, int, double);
	virtual void solve(lexer*,fdm2D*, ghostcell*, vec2D&, vec2D&, int, int, int&, int, double, cpt2D&);
	virtual void setup(lexer*,fdm2D*, ghostcell*,int, cpt2D&);
    
private:
    
    void precon_solve(lexer*,fdm2D*, ghostcell*, vec2D&, vec2D&, int, int, int&, int, double, cpt2D&);
	void precon_setup(lexer*,fdm2D*, ghostcell*,int, cpt2D&);
    

    int num_iterations;
    double final_res_norm;
	int stencil_indices[7];
	int nentries;
   
	int numiter,count,q;
    
    
    //cg
    vec2D sj,rj,r0,vj,tj,pj,precoeff,ph,sh;
	
	int *sizeM,*range;

	const double epsi;

	int margin;
	
	double alpha,beta,w1,w2,w,residual,norm_vj,norm_r0,norm_sj,norm_rj ;
    double r_j1, r_j, sigma ;

};

#endif

#endif

