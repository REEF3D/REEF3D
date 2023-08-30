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

void wave_lib_hdc::time_interpol(lexer *p)
{
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    E[i][j] = E1[i][j]*t1 + E2[i][j]*t2;

    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    U[i][j][k] = U1[i][j][k]*t1 + U2[i][j][k]*t2;
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    V[i][j][k] = V1[i][j][k]*t1 + V2[i][j][k]*t2;
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    W[i][j][k] = W1[i][j][k]*t1 + W2[i][j][k]*t2;
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    Z[i][j][k] = Zsig[k]*(E[i][j]+p->wd-B[i][j]) + B[i][j];
}
