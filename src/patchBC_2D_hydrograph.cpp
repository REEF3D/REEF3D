/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"patchBC_2D.h"
#include"lexer.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC_2D::patchBC_hydrograph_read(lexer *p, ghostcell *pgc, int ID)
{
    char name[100];
	double val;
	int count;
	
    if(ID<10)
	sprintf(name,"hydrograph_Q_000%i.dat",ID);
    
    if(ID<100 && ID>9)
	sprintf(name,"hydrograph_Q_00%i.dat",ID);
    
    if(ID<1000 && ID>99)
	sprintf(name,"hydrograph_Q_0%i.dat",ID);
    
    if(ID<10000 && ID>999)
	sprintf(name,"hydrograph_Q_%i.dat",ID);
    

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
    
    
    patch[qq]->hydroQ_count=count;
	
	p->Darray(patch[qq]->hydroQ,patch[qq]->hydroQ_count,2);
	
	hg.open (name, ios_base::in);
	
	count=0;
	while(!hg.eof())
	{
	hg>>patch[qq]->hydroQ[count][0]>>patch[qq]->hydroQ[count][1];
	++count;
	}
    
}


double patchBC_2D::patchBC_hydrograph_ipol(lexer *p, ghostcell *pgc, int ID)
{

    return 0.0;

}