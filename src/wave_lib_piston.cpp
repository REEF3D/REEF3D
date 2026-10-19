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

#include"wave_lib_piston.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<fstream>

wave_lib_piston::wave_lib_piston(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: piston wavemaker theory"<<endl;
    }
	
    timecount_old=0;
    wtime_old=0.0;
	timecount=1;
	
	setup_done=0;
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
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    double vel;
    
    if(p->wavetime<ts || p->wavetime>te)
	return 0.0;
	
    if((p->wavetime>kinematics[timecount][0]))
    timecount_old=timecount;
    
	while(p->wavetime>kinematics[timecount][0])
	++timecount;

    // restart
    if(p->wavetime-wtime_old > 10.0*p->dt)
    timecount_old = MAX(timecount-1,0);

    wtime_old = p->wavetime;
	
	vel = p->B118*(kinematics[timecount][1]-kinematics[timecount_old][1])/(kinematics[timecount][0]-kinematics[timecount_old][0]);
    
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
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    double vel;
    
    vel = 0.0;

    return vel;
}

double wave_lib_piston::wave_eta(lexer *p, double x, double y)
{
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    double eta=0.0;
    
    double depth = p->wd - p->bed[IJ];
    
    eta =  (kinematics[timecount][1]-kinematics[timecount_old][1])/(kinematics[timecount][0]-kinematics[timecount_old][0]) * sqrt(depth/9.81);    
    
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

void wave_lib_piston::wavemaker_setup(lexer *p, ghostcell *pgc, double depth)
{
    if(setup_done==1)
    return;
    
    setup_done=1;
    
    if(p->B113>0 || p->B119==1)
    if(p->mpirank==0)
    cout<<"Wave_Lib: piston set-up, still-water depth at paddle: "<<depth<<"  B118: "<<p->B118<<endl;
    
    // moving-paddle BC terms from the first-order signal (before correcting it)
    if(p->B119==1 && (p->A10==3 || p->A10==5))
    wm2.make_Qtable(p,pgc,kinematics,ptnum,depth,1,0.0,depth);
    
    // 2nd-order wavemaker correction
    if(p->B113>0)
    wm2.correct(p,pgc,kinematics,ptnum,depth,1,0.0,depth);
}

double wave_lib_piston::wave_paddle_Q(lexer *p, double z)
{
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    if(p->B119!=1 || p->wavetime<ts || p->wavetime>te)
    return 0.0;
    
    return cosgamma*wm2.paddle_Q(p->wavetime,z);
}

// active absorption target (FNPF B 98 4/7): linear surface elevation of the paddle signal at x0 from the paddle,
// tabulated on the first call (and again if x0 changes)
double wave_lib_piston::wave_eta_arc(lexer *p, double x0)
{
    if(setup_done==0 && p->A10!=5)
    wavemaker_setup(p,nullptr,wdt);
    
    if(!wm2.eta_target_ready() || fabs(x0-wm2.eta_target_x0())>1.0e-9)
    wm2.make_etatable(p,nullptr,kinematics,ptnum,wdt,1,0.0,wdt,x0);
    
    if(p->wavetime<ts || p->wavetime>te)
    return 0.0;
    
    return wm2.eta_target(p->wavetime);
}

void wave_lib_piston::wave_prestep(lexer *p, ghostcell *pgc)
{
}
