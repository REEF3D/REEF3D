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

#include"wave_lib_wcp.h"
#include"lexer.h"

void wave_lib_wcp::read_header(lexer *p, ghostcell *pgc)
{
    int is,js,xs,ys;
    
    filename_header(p,pgc);
    
    ifstream header;
	header.open(name, ios::binary);
    
    // read header
    header.open(name, ios::binary);
    
    // origin_ij
        header.read((char*)&iin, sizeof (int));
        is=iin;
        
        header.read((char*)&iin, sizeof (int));
        js=iin;
        
        // origin_xy
        header.read((char*)&ffn, sizeof (float)); 
        xs=ffn;
        
        header.read((char*)&ffn, sizeof (float)); 
        ys=ffn;
        
        // Nx,Ny,Nz
        header.read((char*)&iin, sizeof (int));
        Nx=iin;
        
        header.read((char*)&iin, sizeof (int));
        Ny=iin;
        
        header.read((char*)&iin, sizeof (int));
        Nz=iin;
        
        
        // allocate arrays
        p->Darray(X,Nx);
        p->Darray(Y,Ny);
        p->Darray(Z,Nz);
        p->Darray(B,Nx,Ny);
        
        
        // write coordinates
        for(i=0; i<Nx; ++i)
        {
        header.read((char*)&ffn, sizeof (float)); 
        X[i]=ffn;
        }
            
        for(j=0; j<Ny; ++j)
        {
        header.read((char*)&ffn, sizeof (float)); 
        Y[j]=ffn;
        }
        
        for(k=0; k<Nz; ++k)
        {
        header.read((char*)&ffn, sizeof (float)); 
        Z[k]=ffn;
        }
        
        for(i=0; i<Nx; ++i)
        for(j=0; j<Ny; ++j)
        {
        header.read((char*)&ffn, sizeof (float)); 
        B[i][j]=ffn;
        }  
        
        
        // numer of iterations
        header.read((char*)&iin, sizeof (int));
        numiter=iin;
        
        p->Darray(simtime,numiter);
        
        for(n=0;n<numiter;++n)
        {
        header.read((char*)&iin, sizeof (int));
        simtime[n]=iin;
        }
        
    header.close();
    
}