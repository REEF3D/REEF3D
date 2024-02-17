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

#include"flowfile_out.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void flowfile_out::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    filecount=0;
    
    if(p->mpirank==0)
	mkdir("./REEF3D_FlowFile",0777);
	
	if(p->mpirank==0 && p->P230>0)
	cout<<"FlowFile: "<<probenum<<endl;

	fileout = new ofstream[p->P230];
    headerout = new ofstream[p->P230];

    
    Ni = 1;
    Nj = 1;
    Nk = p->knoz;

    elnum = Ni*Nj*Nk;
    
    
    p->Darray(U,p->P230,elnum);
	p->Darray(V,p->P230,elnum);
	p->Darray(W,p->P230,elnum);
	p->Darray(P,p->P230,elnum);
	p->Darray(K,p->P230,elnum);
	p->Darray(E,p->P230,elnum);
	p->Darray(VT,p->P230,elnum);
	p->Darray(LS,p->P230,elnum);


	p->Iarray(flag,p->P230,elnum);
    p->Iarray(iloc,p->P230);
    
    for(n=0;n<p->P230;++n)
    iloc[n] = p->posf_i(p->P230_x[n]);
    
}

void flowfile_out::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{

	
	for(n=0;n<p->P230;++n)
	for(k=0;k<p->knoz;++k)
	flag[n][k]=0;
    
    i=0;
    j=0;
    for(n=0;n<p->P230;++n)
    {
    i=iloc[n];
   
        for(k=0;k<p->knoz;++k)
        PCHECK
        flag[n][k]=1;
    }
}

