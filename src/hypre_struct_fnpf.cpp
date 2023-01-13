/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"hypre_struct_fnpf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

hypre_struct_fnpf::hypre_struct_fnpf(lexer* p,ghostcell *pgc, int solve_input, int precon_input) : solve_type(solve_input), precon_type(precon_input)
{	
    int vecsize=p->knox*p->knoy*p->knoz; 
    
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(ilower,3);
    p->Iarray(iupper,3);

    if(p->j_dir==1)
    make_grid(p,pgc);	
    
    if(p->j_dir==0)
    make_grid_2Dvert(p,pgc);
}

hypre_struct_fnpf::~hypre_struct_fnpf()
{
}

void hypre_struct_fnpf::startF(lexer* p, ghostcell* pgc, double *f, double *rhs, double *M, int var)
{
    start_solver8(p,pgc,f,rhs,M);
}


void hypre_struct_fnpf::start_solver8(lexer* p, ghostcell* pgc, double *f, double *rhs, double *M)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    if(p->j_dir==1)
    fill_matrix8(p,pgc,f,rhs,M);
    
    if(p->j_dir==0)
    fill_matrix8_2Dvert(p,pgc,f,rhs,M);

    solve(p,pgc);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec8(p,f,rhs,M);
	
	delete_solver5(p,pgc);
}

#endif
