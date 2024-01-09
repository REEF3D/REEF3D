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

#include"flowfile_in.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>


void flowfile_in::header_read(lexer *p, ghostcell *pgc)
{
    // Open File
	if(p->P14==0)
    sprintf(name,"REEF3D-flowheader-%i.r3d",p->I230);
			
    if(p->P14==1)
    sprintf(name,"./REEF3D_FlowFile/REEF3D-flowheader-%i.r3d",p->I230);
	
    
    // count entries
	headerfile.open(name, ios::binary);
    
    entrycount=0;
    
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&ddn, sizeof (double));
    headerfile.read((char*)&iin, sizeof (int));
    headerfile.read((char*)&iin, sizeof (int));
    headerfile.read((char*)&iin, sizeof (int));

    while(!headerfile.eof())
	{
    headerfile.read((char*)&iin, sizeof (int));
    headerfile.read((char*)&ddn, sizeof (double));
    ++entrycount;
    }
    headerfile.close();

    //alloocate arrays
    p->Iarray(iter,entrycount);
    p->Darray(simtime,entrycount);
    
    // count entries
    headerfile.open(name, ios::binary);
    
    q=0;
    headerfile.read((char*)&ddn, sizeof (double));
    xs=ddn;
    headerfile.read((char*)&ddn, sizeof (double));
    xe=ddn;
    headerfile.read((char*)&ddn, sizeof (double));
    ys=ddn;
    headerfile.read((char*)&ddn, sizeof (double));
    ye=ddn;
    headerfile.read((char*)&ddn, sizeof (double));
    zs=ddn;
    headerfile.read((char*)&ddn, sizeof (double));
    ze=ddn;
    headerfile.read((char*)&ddn, sizeof (double));
    deltax=ddn;
    headerfile.read((char*)&iin, sizeof (int));
    Ni=iin;
    headerfile.read((char*)&iin, sizeof (int));
    Nj=iin;
    headerfile.read((char*)&iin, sizeof (int));
    Nk=iin;


//cout<<p->mpirank<<" Ni: "<<Ni<<" Nj: "<<Nj<<" Nk: "<<Nk<<endl;

    while(!headerfile.eof())
	{
    headerfile.read((char*)&iin, sizeof (int));
    iter[q] = iin;
    
    headerfile.read((char*)&ddn, sizeof (double));
    simtime[q] = double(ddn);
    
    //cout<<p->mpirank<<" IT: "<<iter[q]<<" time: "<<simtime[q]<<endl;
    ++q;
    }
    headerfile.close();

    // allocate data
    p->Darray(X0,Nj,Nk);
    p->Darray(Y0,Nj,Nk);
    p->Darray(Z0,Nj,Nk);
    p->Darray(U0,Nj,Nk);
    p->Darray(V0,Nj,Nk);
    p->Darray(W0,Nj,Nk);
    p->Darray(LS0,Nj,Nk);
    p->Darray(P0,Nj,Nk);
    
    
    p->Darray(X1,Nj,Nk);
    p->Darray(Y1,Nj,Nk);
    p->Darray(Z1,Nj,Nk);
    p->Darray(U1,Nj,Nk);
    p->Darray(V1,Nj,Nk);
    p->Darray(W1,Nj,Nk);
    p->Darray(LS1,Nj,Nk);
    p->Darray(P1,Nj,Nk);
    
    
    t0 = simtime[0];
    t1 = simtime[1];
    q0 = iter[0];
    q1 = iter[1];
    
    // Open File 0
	if(p->P14==0)
    sprintf(name0,"REEF3D-flowfile-%i-%i.r3d",p->I230,q0);
			
    if(p->P14==1)
    sprintf(name0,"./REEF3D_FlowFile/REEF3D-flowfile-%i-%i.r3d",p->I230,q0);
    
    flowfile0.open(name0, ios::binary);

    // Open File 1
	if(p->P14==0)
    sprintf(name1,"REEF3D-flowfile-%i-%i.r3d",p->I230,q1);
			
    if(p->P14==1)
    sprintf(name1,"./REEF3D_FlowFile/REEF3D-flowfile-%i-%i.r3d",p->I230,q1);
    
    flowfile1.open(name1, ios::binary);

}
