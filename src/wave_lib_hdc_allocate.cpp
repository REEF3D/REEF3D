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

#include"wave_lib_hdc.h"
#include"lexer.h"

void wave_lib_hdc::allocate(lexer *p, ghostcell *pgc)
{
    p->Darray(U1,Nx,Ny,Nz);
    p->Darray(U2,Nx,Ny,Nz);
    p->Darray(U,Nx,Ny,Nz);
    
    p->Darray(V1,Nx,Ny,Nz);
    p->Darray(V2,Nx,Ny,Nz);
    p->Darray(V,Nx,Ny,Nz);
    
    p->Darray(W1,Nx,Ny,Nz);
    p->Darray(W2,Nx,Ny,Nz);
    p->Darray(W,Nx,Ny,Nz);
    
    p->Darray(E1,Nx,Ny);
    p->Darray(E2,Nx,Ny);
    p->Darray(E,Nx,Ny);
}
