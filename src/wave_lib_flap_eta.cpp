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

#include"wave_lib_flap_eta.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<fstream>

wave_lib_flap_eta::wave_lib_flap_eta(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: flap_eta wavemaker theory; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<endl;
    }
	
	timecount=0;
	
	read(p,pgc);
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));	
}

wave_lib_flap_eta::~wave_lib_flap_eta()
{
}

double wave_lib_flap_eta::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_flap_eta::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_flap_eta::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel,zcoor,fac;
    
    zcoor=p->pos_z();

    if(p->simtime<ts || p->simtime>te || timecount>=ptnum-1)
	return 0.0;
	
	if(p->simtime>eta[timecount+1][0])
	++timecount;
    
    fac = 2.0*(z-p->B111_zs)/(p->B111_ze-p->B111_zs);
	
	vel = sqrt(9.81/wdt) *fac* wave_eta(p,x,y);

    return vel;
}

double wave_lib_flap_eta::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
    
    vel = 0.0;

    return vel;
}

double wave_lib_flap_eta::wave_eta(lexer *p, double x, double y)
{
    double val=0.0;
    
    if(p->simtime<ts || p->simtime>te || timecount>=ptnum-1)
	return 0.0;
    
    val =  ((eta[timecount+1][1]-eta[timecount][1])/(eta[timecount+1][0]-eta[timecount][0]))
            *((p->simtime)-eta[timecount][0]) + eta[timecount][1];
	
    return val;
}

double wave_lib_flap_eta::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_flap_eta::parameters(lexer *p, ghostcell *pgc)
{

}

void wave_lib_flap_eta::read(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1;
	int count;
	
	sprintf(name,"wavemaker_eta.dat");

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	{
		cout<<endl<<("no 'wavemaker_eta.dat' file found")<<endl<<endl;

	}
	
	count=0;
	while(!file.eof())
	{
	file>>val0>>val1;
	if(val0>=p->B117)
	++count;
	}
	
	file.close();

    
    ptnum=count;
	
	p->Darray(eta,ptnum,2);
	
	file.open ("wavemaker_eta.dat", ios_base::in);
	
	count=0;
	while(!file.eof())
	{
	
	file>>val0>>val1;
	
	if(val0>=p->B117)
	{
	eta[count][0] = val0-p->B117;
	eta[count][1] = val1;
	++count;
	}
	}
	
	ts = eta[0][0];
	te = eta[ptnum-1][0];
	
	
	    
}

void wave_lib_flap_eta::wave_prestep(lexer *p, ghostcell *pgc)
{
}
