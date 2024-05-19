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

#include"patchBC.h"
#include"lexer.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_hydrograph_FSF_read(lexer *p, ghostcell *pgc, int qq, int ID)
{
    char name[100];
	double val;
	int count;

    
    sprintf(name,"hydrograph_FSF_%i.dat",ID);
    
// open file------------
	ifstream hg(name, ios_base::in);
	
	if(!hg)
	{
		cout<<endl<<"no "<<name<<" file found !!!"<<endl<<endl;
	}
	
	count=0;
	while(!hg.eof())
	{
	hg>>val;
	++count;
	}
	
	hg.close();
    

	count/=2;
    
    patch[qq]->hydroFSF_count=count;
	
	p->Darray(patch[qq]->hydroFSF, patch[qq]->hydroFSF_count, 2);
	
	hg.open (name, ios_base::in);
	
	count=0;
	while(!hg.eof())
	{
	hg>>patch[qq]->hydroFSF[count][0]>>patch[qq]->hydroFSF[count][1];
	++count;
	}
    
    hg.close();
}

double patchBC::patchBC_hydrograph_FSF_ipol(lexer *p, ghostcell *pgc, int qq, int ID)
{
	double val;
    
    for(int n=0;n<patch[qq]->hydroFSF_count-1;++n)
    if(p->simtime>=patch[qq]->hydroFSF[n][0] && p->simtime<patch[qq]->hydroFSF[n+1][0])
	{
    val = ((patch[qq]->hydroFSF[n+1][1]-patch[qq]->hydroFSF[n][1])/(patch[qq]->hydroFSF[n+1][0]-patch[qq]->hydroFSF[n][0]))*(p->simtime-patch[qq]->hydroFSF[n][0]) + patch[qq]->hydroFSF[n][1];
	}
    
    if(p->count==0 )
    val = patch[qq]->hydroFSF[0][1];
	
	if(p->simtime>=patch[qq]->hydroFSF[patch[qq]->hydroFSF_count-1][0])
	val=patch[qq]->hydroFSF[patch[qq]->hydroFSF_count-1][1];

	return val;
}
