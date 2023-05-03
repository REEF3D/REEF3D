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

void wave_lib_hdc::read_result_continuous(lexer *p, ghostcell *pgc, double **E0, double ***U0, double ***V0, double ***W0, int q0)
{

    // read file_iter
    result.read((char*)&iin, sizeof (int));
    file_iter=iin;
    
    if(p->mpirank==0)
    cout<<"HDC file_iter: "<<file_iter<<endl;
    
    // read
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    {
        result.read((char*)&ffn, sizeof (float));
        E0[i][j]=double(ffn);
    } 
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    {
        result.read((char*)&ffn, sizeof (float)); 
        U0[i][j][k]=double(ffn);
    } 
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    {
        result.read((char*)&ffn, sizeof (float)); 
        V0[i][j][k]=double(ffn);
    } 
    
    for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
    for(k=0; k<Nz; ++k)
    {
        result.read((char*)&ffn, sizeof (float)); 
        W0[i][j][k]=double(ffn);
    }  
}


        
