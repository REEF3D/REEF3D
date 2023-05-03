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

#include"wave_lib_flap_double.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<fstream>

wave_lib_flap_double::wave_lib_flap_double(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: double-hinged flap wavemaker theory; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<endl;
    }
	
	timecount_old=0;
    timecount=1;
    timecount_z=0;
	
	read(p,pgc);
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_flap_double::~wave_lib_flap_double()
{
}

double wave_lib_flap_double::wave_u(lexer *p, double x, double y, double z)
{
   double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_flap_double::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_flap_double::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel,fac,dX;
    
    z+=p->wd;

    if(p->simtime<ts || p->simtime>te)
	return 0.0;
	
	if(z<p->B112_zs || z>p->B112_ze)
	return 0.0;
	
	if((p->simtime>kinematics[timecount][0]))
    timecount_old=timecount;
    
	while(p->simtime>kinematics[timecount][0])
	++timecount;
    
    
	dX=0.0;
    if(z>p->B112_zs && z<p->B112_ze)
    {
    fac = (z-p->B112_zs)/(p->B112_z2-p->B112_zs);
	dX = fac*((kinematics[timecount][1]-kinematics[timecount_old][1])/(kinematics[timecount][0]-kinematics[timecount_old][0]));
    }
    
    if(z>=p->B112_z2 && z<p->B112_ze)
    {
    fac = (z-p->B112_z2)/(p->B112_ze-p->B112_z2);
	dX += fac*((kinematics[timecount][2]-kinematics[timecount_old][2])/(kinematics[timecount][0]-kinematics[timecount_old][0]));
    }
    
    vel = dX;

    return vel;
}

double wave_lib_flap_double::wave_w(lexer *p, double x, double y, double z)
{
    double vel,fac,dZ;
    
    if(p->B115==0)
    return 0.0;
	
	z+=p->wd;

    if(p->simtime<ts || p->simtime>te)
	return 0.0;
	
	if(z<p->B112_zs || z>p->B112_ze)
	return 0.0;
	
	
	dZ=0.0;
    if(z>p->B112_zs && z<p->B112_ze)
    {
    fac = (z-p->B112_zs)/(p->B112_z2-p->B112_zs);
	dZ = fac*((kinematics[timecount][3]-kinematics[timecount_old][3])/(kinematics[timecount][0]-kinematics[timecount_old][0]));
    }
    
    if(z>=p->B112_z2 && z<p->B112_ze)
    {
    fac = (z-p->B112_z2)/(p->B112_ze-p->B112_z2);
	dZ += fac*((kinematics[timecount][4]-kinematics[timecount_old][4])/(kinematics[timecount][0]-kinematics[timecount_old][0]));
    }
    
    vel = dZ;

    return vel;
}

double wave_lib_flap_double::wave_eta(lexer *p, double x, double y)
{
    double eta;

    eta =  0.0;

    return eta;
}

double wave_lib_flap_double::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;

    fi = (x-p->global_xmin)*wave_horzvel(p,x,y,z);
    
    //cout<<"FI: "<<wave_horzvel(p,x,y,z)<<endl;
    
    return fi;
}

void wave_lib_flap_double::parameters(lexer *p, ghostcell *pgc)
{

}

void wave_lib_flap_double::read(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1,val2;
    double sign1,sign2;
    double beta,s,sign;
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
	file>>val0>>val1>>val2;
	if(val0>=p->B117)
	++count;
	}
	
	file.close();

    
    ptnum=count;
	
	p->Darray(kinematics,ptnum,5);
	
	file.open ("wavemaker.dat", ios_base::in);
	
	count=0;
	while(!file.eof())
	{
	
        file>>val0>>val1>>val2;
        
        if(val0>=p->B117)
        {
        kinematics[count][0] = val0-p->B117;
        kinematics[count][1] = val1;
        kinematics[count][2] = val2;
        ++count;
        }
	}
	
	ts = kinematics[0][0];
	te = kinematics[ptnum-1][0];
    
    // convert from angles to X
    if(p->B116==2)
    for(int qn=0; qn<ptnum; ++qn)
    {
    sign1 = -(kinematics[qn][1]/(fabs(kinematics[qn][1])>1.0e-20?fabs(kinematics[qn][1]):1.0e20));
    sign2 = -(kinematics[qn][2]/(fabs(kinematics[qn][2])>1.0e-20?fabs(kinematics[qn][2]):1.0e20));
    
    kinematics[qn][1] = sign1*fabs(sin(kinematics[qn][1])*(p->B112_z2-p->B112_zs));
    kinematics[qn][2] = sign2*fabs(sin(kinematics[qn][2])*(p->B112_ze-p->B112_z2));
    }
    
    // calculate vertical component
    for(int qn=0; qn<ptnum; ++qn)
    {
    sign = -1.0;
    
    beta = asin(kinematics[qn][1]/(p->B112_z2-p->B112_zs));
    s = 2.0*(p->B112_z2-p->B112_zs)*sin(0.5*beta);
    kinematics[qn][3] = sign*sqrt(pow(s,2.0) - pow(kinematics[qn][1],2.0));
    
    beta = asin(kinematics[qn][2]/(p->B112_ze-p->B112_z2));
    s = 2.0*(p->B112_ze-p->B112_z2)*sin(0.5*beta);
    kinematics[qn][4] = sign*sqrt(pow(s,2.0) - pow(kinematics[qn][2],2.0));
    }
    
    /*
    cout<<"flap kinematics"<<endl;
    if(p->mpirank==0)
    for(int qn=0; qn<ptnum; ++qn)
    cout<<kinematics[qn][1]<<" "<<kinematics[qn][2]<<" | "<<kinematics[qn][3]<<" "<<kinematics[qn][4]<<endl;
    */
	    
}

void wave_lib_flap_double::wave_prestep(lexer *p, ghostcell *pgc)
{
}
