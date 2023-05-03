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

#include"potentialfile_in.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void potentialfile_in::read0(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{    
    
    
    // Open File 
	if(p->P14==0)
    sprintf(name0,"REEF3D-potentialfile-%i-%i.r3d",p->I230,q0);
			
    if(p->P14==1)
    sprintf(name,"./REEF3D_PotentialFile/REEF3D-flowfile-%i-%i.r3d",p->I230,q0);
    
    
    // read header part
    potentialfile.open(name, ios::binary);

    while(!potentialfile.eof())
	{
    potentialfile.read((char*)&ffn, sizeof (float));
    bedlevel = ffn;
    
    potentialfile.read((char*)&ffn, sizeof (float));
    waterlevel = ffn;
    
    potentialfile.read((char*)&iin, sizeof (int));
    Nk = iin;
    
    p->Darray(S,Nk);
    
        for(qn=0; qn<Nk; ++qn)
        {
        potentialfile.read((char*)&ffn, sizeof (float));
        S[qn]=ffn;
        }
    }
    potentialfile.close();
    
    
    
    
    
    
    // count entries
    potentialfile.open(name, ios::binary);
    
    entrycount=0;
    
    potentialfile.read((char*)&ffn, sizeof (float));
    potentialfile.read((char*)&ffn, sizeof (float));
    potentialfile.read((char*)&iin, sizeof (int));
    
    for(qn=0; qn<Nk; ++qn)
    potentialfile.read((char*)&ffn, sizeof (float));

    
    while(!potentialfile.eof())
	{
    for(qn=0; qn<2+2*Nk; ++qn)
    potentialfile.read((char*)&ffn, sizeof (float));
    
    ++entrycount;
    }
    potentialfile.close();
    
    
    
    // allocate data
    p->Darray(T,entrycount);
    p->Darray(E,entrycount);
    
    p->Darray(U,entrycount,Nk);
    p->Darray(W,entrycount,Nk);
    
    
    t0 = T[0];
    t1 = T[1];

    
    // start reading 
    potentialfile.read((char*)&ffn, sizeof (float));
    xloc[n]=ffn;



/*
        ffn = float(c->bed(i,j));
        fileout[n].write((char*)&ffn, sizeof (float));
        
        iin=p->knoz;
        fileout[n].write((char*)&iin, sizeof (int));
        
        for(qn=0;qn<p->knoz+1;++qn)
        {
        ffn = float(p->ZN[KP]);
        fileout[n].write((char*)&ffn, sizeof (float));
        }
*/
    

    
   


 
    /*
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
    }*/
}
    
