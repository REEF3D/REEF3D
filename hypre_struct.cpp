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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

hypre_struct::hypre_struct(lexer* p,fdm* a,ghostcell *pgc) : cval4(p)
{	
    int vecsize=p->knox*p->knoy*p->knoz; 
    
    if(p->A10==3 && p->A300==1)
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(ilower,3);
    p->Iarray(iupper,3);
    p->Darray(values,vecsize*7);
    
    if(p->j_dir==1)
    make_grid(p,a,pgc);	
    
    if(p->j_dir==0)
    make_grid_2Dvert(p,a,pgc);
    
    
    count=0;
    FLUIDLOOP
    {
    cval4(i,j,k)=count;
    ++count;
    }
}

hypre_struct::~hypre_struct()
{
}

void hypre_struct::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
    if(var>=1 && var<=4)
    start_solver1234(p,a,pgc,f,xvec,rhsvec,var);
    
    if(var==5)
    start_solver5(p,a,pgc,f,xvec,rhsvec,var);
}

void hypre_struct::startF(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var, int gcv, double stop_crit)
{
    if(var==7)
    start_solver7(p,pgc,f,rhs,M,var);
    
    if(var==8)
    start_solver8(p,c,pgc,f,rhs,M,var);
}

void hypre_struct::start_solver1234(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var)
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
    fillbackvec1(p,f,xvec,var);
    
    if(var==2)
    fillbackvec2(p,f,xvec,var);
    
    if(var==3)
    fillbackvec3(p,f,xvec,var);
    
    if(var==4)
    fillbackvec4(p,f,xvec,var);
	
	delete_solver1234(p,pgc);
}

void hypre_struct::start_solver5(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var)
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
        
    fillbackvec4(p,f,xvec,var);
	
	delete_solver5(p,pgc);
}


void hypre_struct::start_solver7(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix7(p,pgc,f,rhs,M);

    solve(p,pgc);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec7(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::start_solver8(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
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
        
    fillbackvec8(p,c,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void hypre_struct::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
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
