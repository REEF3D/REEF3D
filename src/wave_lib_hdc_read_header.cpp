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

#include"wave_lib_hdc.h"
#include"lexer.h"

void wave_lib_hdc::read_header(lexer *p, ghostcell *pgc)
{
    int is,js,xs,ys;
    ifstream header;
    
    // filename
    filename_header(p,pgc);
    
    // open header
    header.open(name, ios::binary);
    
        // file_version
        header.read((char*)&iin, sizeof (int));
        file_version=iin;
        
        // file_type
        header.read((char*)&iin, sizeof (int));
        file_type=iin;
        
        // jdir
        header.read((char*)&iin, sizeof (int));
        jdir=iin;
        
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
        p->Darray(Zsig,Nz);
        p->Darray(Z,Nx,Ny,Nz);
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
        Zsig[k]=ffn;
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
        
        // delta iterations
        header.read((char*)&iin, sizeof (int));
        diter=iin;
        
        
        p->Darray(simtime,numiter);
        
        for(n=0;n<numiter;++n)
        {
        header.read((char*)&ddn, sizeof (double));
        simtime[n]=ddn;
        }
        
        t_start = simtime[0];
        t_end   = simtime[numiter-1];

        
        if(Nx==0)
        {
        Xstart=1.0;
        Xend=-1.0;
        }
        
        if(Nx>0)
        {
        Xstart = X[0];
        Xend = X[Nx-1];
        }
        
        
        
        if(Ny==0)
        {
        Ystart=1.0;
        Yend=-1.0;
        }
        
        if(Ny>0)
        {
        Ystart = Y[0];
        Yend = Y[Ny-1];
        }
        
        
    header.close();
    
}
