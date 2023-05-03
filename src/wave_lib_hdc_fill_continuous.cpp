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

void wave_lib_hdc::fill_result_continuous(lexer *p, ghostcell *pgc)
{

    // fill
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    E1[i][j]=E2[i][j];
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    U1[i][j][k]=U2[i][j][k];
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    V1[i][j][k]=V2[i][j][k];

    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    W1[i][j][k]=W2[i][j][k];
    
    
}


        
