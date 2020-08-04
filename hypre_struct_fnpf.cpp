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
--------------------------------------------------------------------*/

#include"hypre_struct_fnpf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

hypre_struct_fnpf::hypre_struct_fnpf(lexer* p,fdm* a,ghostcell *pgc, int solve_input, int precon_input) : solve_type(solve_input), precon_type(precon_input)
{	
    int vecsize=p->knox*p->knoy*p->knoz; 
    
    if(p->A10==3)
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(ilower,3);
    p->Iarray(iupper,3);

    if(p->j_dir==1)
    make_grid(p,a,pgc);	
    
    if(p->j_dir==0)
    make_grid_2Dvert(p,a,pgc);
}

hypre_struct_fnpf::~hypre_struct_fnpf()
{
}

void hypre_struct_fnpf::startF(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, double *rhs, double *M, int var, int gcv, double stop_crit)
{
    start_solver8(p,c,pgc,f,rhs,M);
}


void hypre_struct_fnpf::start_solver8(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, double *rhs, double *M)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    if(p->j_dir==1)
    fill_matrix8(p,c,pgc,f,rhs,M);
    
    if(p->j_dir==0)
    fill_matrix8_2Dvert(p,c,pgc,f,rhs,M);

    solve(p,pgc);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec8(p,c,f,rhs,M);
	
	delete_solver5(p,pgc);
}

#endif
