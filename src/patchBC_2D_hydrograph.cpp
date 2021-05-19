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
   /* char name[100];
	double val;
	int count;
	
    
	sprintf(name,"hydrograph.dat");
    
    if(num<10)
	sprintf(name,"REEF3D-CFD-00000%i.pvtu",num);

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
    
    
    hydro_in_count=count;
	
	p->Darray(hydro,hydro_count,2);
	
	hg.open ("hydrograph.dat", ios_base::in);
	
	count=0;
	while(!hg.eof())
	{
	hg>>hydro_in[count][0]>>hydro_in[count][1];
	++count;
	}
    */
}


double patchBC_2D::patchBC_hydrograph_ipol(lexer *p, ghostcell *pgc, int ID)
{

    return 0.0;

}