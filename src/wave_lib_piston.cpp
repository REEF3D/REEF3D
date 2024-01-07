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

#include"wave_lib_piston.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<fstream>

wave_lib_piston::wave_lib_piston(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: piston wavemaker theory";
    }
	
    timecount_old=0;
	timecount=1;
	
	read(p,pgc);
	
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
    
}

wave_lib_piston::~wave_lib_piston()
{
}

double wave_lib_piston::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_piston::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_piston::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;
    

    if(p->simtime<ts || p->simtime>te)
	return 0.0;
	
    if((p->simtime>kinematics[timecount][0]))
    timecount_old=timecount;
    
	while(p->simtime>kinematics[timecount][0])
	++timecount;
	
	vel = (kinematics[timecount][1]-kinematics[timecount_old][1])/(kinematics[timecount][0]-kinematics[timecount_old][0]);
    
    if(p->B110==1)
    {
    z+=p->wd;
    
    if(z<p->B110_zs || z>p->B110_ze)
    vel=0.0;
    }
    
    return vel;
}

double wave_lib_piston::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
    
    vel = 0.0;

    return vel;
}

double wave_lib_piston::wave_eta(lexer *p, double x, double y)
{
    double eta;

    eta =  0.0;

    return eta;
}

double wave_lib_piston::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;

    fi = (x-p->global_xmin)*wave_horzvel(p,x,y,z);
    
    return fi;
}

void wave_lib_piston::parameters(lexer *p, ghostcell *pgc)
{

}

void wave_lib_piston::read(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1;
	int count;
	
	sprintf(name,"wavemaker.dat");

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	{
		cout<<endl<<("no 'wavemaker.dat' file found")<<endl<<endl;

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
	
	p->Darray(kinematics,ptnum,2);
	
	file.open ("wavemaker.dat", ios_base::in);
	
	count=0;
	while(!file.eof())
	{
	
	file>>val0>>val1;
	
	if(val0>=p->B117)
	{
	kinematics[count][0] = val0-p->B117;
	kinematics[count][1] = val1;
	++count;
	}
	}
	
	ts = kinematics[0][0];
	te = kinematics[ptnum-1][0];
	    
}

void wave_lib_piston::wave_prestep(lexer *p, ghostcell *pgc)
{
}
