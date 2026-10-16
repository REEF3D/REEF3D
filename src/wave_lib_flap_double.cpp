/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
    cout<<"Wave_Lib: double-hinged flap wavemaker theory"<<endl;
    }
	
	timecount_old=0;
	wtime_old=0.0;
    timecount=1;
    timecount_z=0;
	
	setup_done=0;
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
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    double vel,fac,dX;
    
    z+=p->wd;

    if(p->wavetime<ts || p->wavetime>te)
	return 0.0;
	
	if(z<p->B112_zs || z>p->B112_ze)
	return 0.0;
	
	if((p->wavetime>kinematics[timecount][0]))
    timecount_old=timecount;
    
	while(p->wavetime>kinematics[timecount][0])
	++timecount;

    // restart
    if(p->wavetime-wtime_old > 10.0*p->dt)
    timecount_old = MAX(timecount-1,0);

    wtime_old = p->wavetime;
    
    
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
    
    vel = dX*p->B118;

    return vel;
}

double wave_lib_flap_double::wave_w(lexer *p, double x, double y, double z)
{
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    double vel,fac,dZ;
    
    if(p->B115==0)
    return 0.0;
	
	z+=p->wd;

    if(p->wavetime<ts || p->wavetime>te)
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
    
    vel = dZ*p->B118;

    return vel;
}

double wave_lib_flap_double::wave_eta(lexer *p, double x, double y)
{
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
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
    
    /*
    cout<<"flap kinematics"<<endl;
    if(p->mpirank==0)
    for(int qn=0; qn<ptnum; ++qn)
    cout<<kinematics[qn][1]<<" "<<kinematics[qn][2]<<" | "<<kinematics[qn][3]<<" "<<kinematics[qn][4]<<endl;
    */
	    
}

void wave_lib_flap_double::wavemaker_setup(lexer *p, ghostcell *pgc, double depth)
{
    double beta,s,sign;
    
    if(setup_done==1)
    return;
    
    setup_done=1;
    
    // hinge elevations above the bed at the paddle (B112 is measured with the still water level at wd)
    const double zb = p->wd - depth;
    const double zs = p->B112_zs - zb;
    const double z2 = p->B112_z2 - zb;
    const double ze = p->B112_ze - zb;
    
    if(p->B113>0 || p->B119==1)
    if(p->mpirank==0)
    cout<<"Wave_Lib: double flap set-up, still-water depth at paddle: "<<depth<<"  hinges above bed: "<<zs<<" "<<z2<<"  top above bed: "<<ze<<"  B118: "<<p->B118<<endl;
    
    if((p->B113>0 || p->B119==1) && zs<0.0 && p->mpirank==0)
    cout<<"Wave_Lib: WARNING lower hinge below the bed at the paddle, check B112 and the depth (B94)"<<endl;
    
    // moving-paddle BC terms from the first-order signal (before correcting it)
    if(p->B119==1 && (p->A10==3 || p->A10==5))
    wm2.make_Qtable_double(p,pgc,kinematics,ptnum,depth,zs,z2,ze);
    
    // 2nd-order wavemaker correction, split over the two flaps in proportion to the first-order motion
    if(p->B113>0)
    wm2.correct_double(p,pgc,kinematics,ptnum,depth,zs,z2,ze);
    
    // calculate vertical component
    for(int qn=0; qn<ptnum; ++qn)
    {
    sign = -1.0;
    
    beta = asin(MAX(-1.0,MIN(1.0,kinematics[qn][1]/(p->B112_z2-p->B112_zs))));
    s = 2.0*(p->B112_z2-p->B112_zs)*sin(0.5*beta);
    kinematics[qn][3] = sign*sqrt(MAX(0.0,pow(s,2.0) - pow(kinematics[qn][1],2.0)));
    
    beta = asin(MAX(-1.0,MIN(1.0,kinematics[qn][2]/(p->B112_ze-p->B112_z2))));
    s = 2.0*(p->B112_ze-p->B112_z2)*sin(0.5*beta);
    kinematics[qn][4] = sign*sqrt(MAX(0.0,pow(s,2.0) - pow(kinematics[qn][2],2.0)));
    }
}

double wave_lib_flap_double::wave_paddle_Q(lexer *p, double z)
{
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    if(p->B119!=1 || p->wavetime<ts || p->wavetime>te)
    return 0.0;
    
    return cosgamma*wm2.paddle_Q(p->wavetime,z);
}

void wave_lib_flap_double::wave_prestep(lexer *p, ghostcell *pgc)
{
}
