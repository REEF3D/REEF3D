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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<fstream>


double ioflow_f::hydrograph_ipol(lexer *p, ghostcell* pgc, double ** hydro, int hydrocount)
{
	double val;
    
    for(n=0;n<hydrocount-1;++n)
    if(p->simtime>=hydro[n][0] && p->simtime<hydro[n+1][0])
	{
    val = ((hydro[n+1][1]-hydro[n][1])/(hydro[n+1][0]-hydro[n][0]))*(p->simtime-hydro[n][0]) + hydro[n][1];
	}
    
    if(p->count==0 )
    val = hydro[0][1];
	
	if(p->simtime>=hydro[hydrocount-1][0])
	val=hydro[hydrocount-1][1];
	
	return val;
	
}

void ioflow_f::hydrograph_in_read(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val;
	int count;
	
	sprintf(name,"hydrograph.dat");

// open file------------
	ifstream hg(name, ios_base::in);
	
	if(!hg)
	{
		cout<<endl<<("no 'hydrograph.dat' file found")<<endl<<endl;

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
	
	p->Darray(hydro_in,hydro_in_count,2);
	
	hg.open ("hydrograph.dat", ios_base::in);
	
	count=0;
	while(!hg.eof())
	{
	hg>>hydro_in[count][0]>>hydro_in[count][1];
	++count;
	}
    
}

void ioflow_f::hydrograph_out_read(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val;
	int count;
	
	sprintf(name,"hydrograph_out.dat");

// open file------------
	ifstream hg(name, ios_base::in);
	
	if(!hg)
	{
		cout<<endl<<("no 'hydrograph_out.dat' file found")<<endl<<endl;

	}
	
	count=0;
	while(!hg.eof())
	{
	hg>>val;
	++count;
	}
	
	hg.close();
	
	count/=2;
    
    
    hydro_out_count=count;
	
	p->Darray(hydro_out,hydro_out_count,2);
	
	hg.open ("hydrograph_out.dat", ios_base::in);
	
	count=0;
	while(!hg.eof())
	{
	hg>>hydro_out[count][0]>>hydro_out[count][1];
	++count;
	}
    
}
