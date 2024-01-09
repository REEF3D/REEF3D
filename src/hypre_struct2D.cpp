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

#include"hypre_struct2D.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"ghostcell.h"
#include"field.h"
#include"vec2D.h"

hypre_struct2D::hypre_struct2D(lexer* p,ghostcell *pgc)
{	
    int vecsize=p->knox*p->knoy; 
    
    p->Iarray(ilower,2);
    p->Iarray(iupper,2);
    p->Darray(values,vecsize*5);
    
    make_grid(p,pgc);	  
}

hypre_struct2D::~hypre_struct2D()
{
}

void hypre_struct2D::start(lexer* p, ghostcell* pgc, slice &f, matrix2D &M, vec2D& xvec, vec2D& rhsvec, int var)
{    
	create_solvers(p,pgc);
    
    // fill for cfd
    fill_matrix(p,pgc,M,f,rhsvec);
    
    solve(p,pgc);
        
    fillbackvec(p,f,xvec,var);
	
	delete_solvers(p,pgc);
}


#endif
