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

#include"hypre_sstruct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

hypre_sstruct::hypre_sstruct(lexer* p,fdm* a,ghostcell *pgc)
{	
    int vecsize=p->knox*p->knoy*p->knoz; 
    
    if(p->A10==3 && p->A300==1)
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(ilower,3);
    p->Iarray(iupper,3);
    p->Darray(values,vecsize*13);
    
    if(p->A320!=2)
    make_grid_7p(p,a,pgc);	
    
    if(p->A320==2)
    make_grid_13p(p,a,pgc);	
}

hypre_sstruct::~hypre_sstruct()
{
}

void hypre_sstruct::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
    if(var>=1 && var<=4)
    start_solver1234(p,a,pgc,f,xvec,rhsvec,var);
    
    if(var==5)
    start_solver5(p,a,pgc,f,xvec,rhsvec,var);
}

void hypre_sstruct::startF(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var, int gcv, double stop_crit)
{
    if(var==7)
    start_solver7(p,pgc,f,rhs,M,var);
    
    if(var==8)
    start_solver8(p,pgc,f,rhs,M,var);
    
    if(var==10)
    start_solver10(p,pgc,f,rhs,M,var);
}

void hypre_sstruct::start_solver1234(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var)
{
    numiter=0;
	p->solveriter=0;
	
	create_solver1234(p,pgc);
    
    if(var==1)
    fill_matrix1(p,a,pgc,f);
    
    if(var==2)
    fill_matrix2(p,a,pgc,f);
    
    if(var==3)
    fill_matrix3(p,a,pgc,f);
    
    if(var==4)
    fill_matrix4(p,a,pgc,f);
    
    
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

void hypre_sstruct::start_solver5(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var)
{
	numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix4(p,a,pgc,f);

    solve(p);
	
	p->solveriter=num_iterations;
        
    fillbackvec4(p,f,xvec,var);
	
	delete_solver5(p,pgc);
}


void hypre_sstruct::start_solver7(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix7(p,pgc,f,rhs,M);

    solve(p);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec7(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_sstruct::start_solver8(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix8(p,pgc,f,rhs,M);

    solve(p);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec8(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_sstruct::start_solver10(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix10(p,pgc,f,rhs,M);


    solve(p);

	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec10(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_sstruct::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void hypre_sstruct::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
}

void hypre_sstruct::fillxvec1(lexer* p, fdm* a, field& f)
{
}
	
void hypre_sstruct::fillxvec2(lexer* p, fdm* a, field& f)
{
}
	
void hypre_sstruct::fillxvec3(lexer* p, fdm* a, field& f)
{
}
	
void hypre_sstruct::fillxvec4(lexer* p, fdm* a, field& f)
{
}

#endif
