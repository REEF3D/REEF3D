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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

hypre_struct::hypre_struct(lexer* p,ghostcell *pgc, int solve_input, int precon_input) : solve_type(solve_input), precon_type(precon_input)
{	
    int vecsize=p->knox*p->knoy*p->knoz; 
    
    p->Iarray(CVAL4,p->imax*p->jmax*(p->kmax+2));
    
    if(p->A10==3)
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(ilower,3);
    p->Iarray(iupper,3);
    
    if(p->D33==0)
    p->Darray(values,vecsize*7);
    
    if(p->D33==1)
    p->Darray(values,vecsize*15);
    
    if(p->j_dir==1 && p->D33==0)
    make_grid(p,pgc);	
    
    if(p->j_dir==0 && p->D33==0)
    make_grid_2Dvert(p,pgc);
    
    if(p->j_dir==1 && p->D33==1)
    make_grid_15pt(p,pgc);	
    
    if(p->j_dir==0 && p->D33==1)
    make_grid_2D_9pt(p,pgc);
    
    
    count=0;
    FLUIDLOOP
    {
    CVAL4[IJK]=count;
    ++count;
    }
}

hypre_struct::~hypre_struct()
{
}

void hypre_struct::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& rhsvec, int var)
{
    if(var>=1 && var<=4)
    start_solver1234(p,a,pgc,f,rhsvec,var);

    if(var==5)
    start_solver5(p,a,pgc,f,rhsvec,var);
}

void hypre_struct::startf(lexer* p, ghostcell* pgc, field &f, vec& rhs, matrix_diag &M, int var)
{
    start_solver4f(p, pgc,f,rhs,M,var);
}

void hypre_struct::startF(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    if(var==7)
    start_solver7(p,pgc,f,rhs,M,var);
    
    if(var==8)
    start_solver8(p,pgc,f,rhs,M,var);
    
    if(var==9)
    start_solver9(p,pgc,f,rhs,M,var);
}

void hypre_struct::startV(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    if(var==4)
    start_solver4V(p,pgc,f,rhs,M,var); 
}

void hypre_struct::startM(lexer* p, ghostcell* pgc, double *x, double *rhs, double *M, int var)
{
}

void hypre_struct::start_solver1234(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& rhsvec, int var)
{
    numiter=0;
	p->solveriter=0;
	
	create_solver1234(p,pgc);
    
    if(var==1)
    {
        if(p->j_dir==1)
        fill_matrix1(p,a,pgc,f);
        
        if(p->j_dir==0)
        fill_matrix1_2Dvert(p,a,pgc,f);
    }
    
    if(var==2)
    {
        if(p->j_dir==1)
        fill_matrix2(p,a,pgc,f);
        
        if(p->j_dir==0)
        fill_matrix2_2Dvert(p,a,pgc,f);
    }
    
    if(var==3)
    {
        if(p->j_dir==1)
        fill_matrix3(p,a,pgc,f);
        
        if(p->j_dir==0)
        fill_matrix3_2Dvert(p,a,pgc,f);
    }
    
    if(var==4)
    {
        if(p->j_dir==1)
        fill_matrix4(p,a,pgc,f);
        
        if(p->j_dir==0)
        fill_matrix4_2Dvert(p,a,pgc,f);
    }
    
    
    solve1234(p);
        
    
    if(var==1)
    fillbackvec1(p,f,var);
    
    if(var==2)
    fillbackvec2(p,f,var);
    
    if(var==3)
    fillbackvec3(p,f,var);
    
    if(var==4)
    fillbackvec4(p,f,var);
	
	delete_solver1234(p,pgc);
}

void hypre_struct::start_solver5(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& rhsvec, int var)
{
	numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);
    
    if(p->j_dir==1)
    fill_matrix4(p,a,pgc,f);
    
    if(p->j_dir==0)
    fill_matrix4_2Dvert(p,a,pgc,f);

    solve(p,pgc);
	
	p->solveriter=num_iterations;
        
    fillbackvec4(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::start_solver4V(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);
    
    if(p->j_dir==1)
    fill_matrix4V(p,pgc,f,rhs,M);
    
    if(p->j_dir==0)
    fill_matrix4V_2D(p,pgc,f,rhs,M);

    solve(p,pgc);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec4V(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::start_solver4f(lexer* p, ghostcell* pgc, field &f, vec& rhs, matrix_diag &M, int var)
{
	numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);
    
    if(p->j_dir==1)
    fill_matrix4f(p,pgc,f,rhs,M);
    
    if(p->j_dir==0)
    fill_matrix4f_2Dvert(p,pgc,f,rhs,M);

    solve(p,pgc);
	
	p->solveriter=num_iterations;
        
    fillbackvec4(p,f,var);
	
	delete_solver5(p,pgc);
}


void hypre_struct::start_solver7(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);
    
    if(p->j_dir==1)
    fill_matrix7(p,pgc,f,rhs,M);
    
    if(p->j_dir==0)
    fill_matrix7_2Dvert(p,pgc,f,rhs,M);

    solve(p,pgc);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec7(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::start_solver8(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
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
        
    fillbackvec8(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::start_solver9(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    if(p->j_dir==1)
    fill_matrix9(p,pgc,f,rhs,M);
    
    if(p->j_dir==0)
    fill_matrix9_2Dvert(p,pgc,f,rhs,M);

    solve(p,pgc);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec9(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::fillxvec1(lexer* p, fdm* a, field& f)
{
}
	
void hypre_struct::fillxvec2(lexer* p, fdm* a, field& f)
{
}
	
void hypre_struct::fillxvec3(lexer* p, fdm* a, field& f)
{
}
	
void hypre_struct::fillxvec4(lexer* p, fdm* a, field& f)
{
}

#endif
