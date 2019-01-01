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

#include"flowfile_in.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>
/*
void flowfile_in::read0(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{    
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    X0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    Y0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    Z0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    U0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    V0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    W0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    P0[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ddn, sizeof (double));
    LS0[j][k]=ddn;
    }
}
    
 void flowfile_in::read1(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{      
    // read ff1
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    X1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    Y1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    Z1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    U1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    V1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    W1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    P1[j][k]=ddn;
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ddn, sizeof (double));
    LS1[j][k]=ddn;
    }

}

*/

void flowfile_in::read0(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{    
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    X0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    Y0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    Z0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    U0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    V0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    W0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    P0[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile0.read((char*)&ffn, sizeof (float));
    LS0[j][k]=double(ffn);
    }
}
    
 void flowfile_in::read1(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{      
    // read ff1
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    X1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    Y1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    Z1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    U1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    V1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    W1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    P1[j][k]=double(ffn);
    }
    
    for(j=0;j<Nj;++j)
    for(k=0;k<Nk;++k)
    {
    flowfile1.read((char*)&ffn, sizeof (float));
    LS1[j][k]=double(ffn);
    }

}







