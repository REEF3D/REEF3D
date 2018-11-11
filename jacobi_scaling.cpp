/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"jacobi_scaling.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"cpt.h"

jacobi_scaling::jacobi_scaling(lexer* p,fdm *a, ghostcell *pgc):aii(p),epsi(1.0e-19)
{	
}

jacobi_scaling::~jacobi_scaling()
{
}

void jacobi_scaling::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
	if(var==1)
	sizeM=p->sizeM1;
	
	if(var==2)
	sizeM=p->sizeM2;
	
	if(var==3)
	sizeM=p->sizeM3;
	
	if(var==4)
	sizeM=p->sizeM4;
	
	NLOOP
	aii.V[n]=-1.0/(a->M.p[n]+epsi);
}

void jacobi_scaling::start(lexer* p,fdm* a, ghostcell* pgc, field &xfield, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
}

void jacobi_scaling::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{                    
}
	
void jacobi_scaling::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
	if(var==1)
	NLOOP
	xvec.V[n]=rhsvec.V[n]*aii.V[n];	
	
	if(var==2)
	NLOOP
	xvec.V[n]=rhsvec.V[n]*aii.V[n];
	
	if(var==3)
	NLOOP
	xvec.V[n]=rhsvec.V[n]*aii.V[n];	
	
	if(var==4)
	NLOOP
	xvec.V[n]=rhsvec.V[n]*aii.V[n];
	
	solveriter=1;
}

void jacobi_scaling::fillxvec1(lexer* p, fdm* a, field& f)
{
}

void jacobi_scaling::fillxvec2(lexer* p, fdm* a, field& f)
{
}

void jacobi_scaling::fillxvec3(lexer* p, fdm* a, field& f)
{
}

void jacobi_scaling::fillxvec4(lexer* p, fdm* a, field& f)
{
}



