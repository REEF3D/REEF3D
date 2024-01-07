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

#include"wave_lib_solitary_3rd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_solitary_3rd::wave_lib_solitary_3rd(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: 3rd-order solitary wave "<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_solitary_3rd::~wave_lib_solitary_3rd()
{
}

double wave_lib_solitary_3rd::wave_u(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);
	
    return cosgamma*vel;
}

double wave_lib_solitary_3rd::wave_v(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);
	
    return singamma*vel;
}

double wave_lib_solitary_3rd::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel,eta,e,s,t,aval;
	
	teta = -(wC*(p->simtime) - (x - X0));
	
	e = wH/wdt;
	
	aval = (1.0/wdt)*sqrt(0.75*e) * (1.0 - (5.0/8.0)*e + (71.0/128.0)*e*e);
	
	s = 1.0/cosh(aval*teta);
	
	t = tanh(aval*teta);
		
	vel = sqrt(9.81*wdt)*(e*s*s + e*e*(-0.75*s*s + s*s*t*t + pow((wdt+z)/wdt,2.0)*(0.75*s*s - (9.0/4.0)*s*s*t*t)
						+ e*e*e*((21.0/40.0)*s*s - s*s*t*t - (6.0/5.0)*s*s*s*s*t*t + pow((wdt+z)/wdt,2.0)*(-(9.0/4.0)*s*s + (15.0/4.0)*s*s*t*t + (15.0/2.0)*s*s*s*s*t*t)
								 + pow((wdt+z)/wdt,4.0) * (3.0/8.0)*s*s - (45.0/16.0)*s*s*s*s*t*t)));
		
    return vel;
}

double wave_lib_solitary_3rd::wave_w(lexer *p, double x, double y, double z)
{
    double vel,eta,e,s,t,aval;
	
	teta = -(wC*(p->simtime) - (x - X0));
	
	e = wH/wdt;
	
	aval = (1.0/wdt)*sqrt(0.75*e) * (1.0 - (5.0/8.0)*e + (71.0/128.0)*e*e);
	
	s = 1.0/cosh(aval*teta);
	
	t = tanh(aval*teta);
			
	vel = sqrt(9.81*wdt)*sqrt(3.0*e)*((wdt+z)/wdt)*s*s*t*( e
		+ e*e*((3.0/8.0)*s*s + 2.0*s*s + pow((wdt+z)/wdt,2.0)*((-1.0/2.0) + 1.5*s*s))
		+ e*e*e*((49.0/640.0) -(17.0/20.0)*s*s - (18.0/5.0)*s*s*s*s 
			+ pow((wdt+z)/wdt,2.0)*(-(13.0/16.0) - (25.0/16.0)*s*s + (15.0/2.0)*s*s*s*s)
			+ pow((wdt+z)/wdt,4.0)*(-(3.0/40.0) + (9.0/8.0)*s*s - (27.0/16.0)*s*s*s*s)));
				
    return vel;
}

double wave_lib_solitary_3rd::wave_eta(lexer *p, double x, double y)
{
    double eta,e,t,s,aval;

	
	teta = -(wC*(p->simtime) - (x - X0));
	
	e = wH/wdt;
	
	aval = (1.0/wdt)*sqrt(0.75*e) * (1.0 - (5.0/8.0)*e + (71.0/128.0)*e*e);
	
	s = 1.0/cosh(aval*teta);
	
	t = tanh(aval*teta);
	
	eta =  wdt*(e*s*s - 0.75*e*e*s*s*t*t + e*e*e*((5.0/8.0)*s*s*t*t - (101.0/80.0)*s*s*s*s*t*t));
	
	return eta;	
}

double wave_lib_solitary_3rd::wave_fi(lexer *p, double x, double y, double z)
{
    double fi=0.0;
    
    return fi;
}

void wave_lib_solitary_3rd::parameters(lexer *p, ghostcell *pgc)
{
	double e = wH/wdt;
	
	wC = sqrt(9.81*wdt)*sqrt(1.0 + e  - (1.0/20.0)*e*e - (3.0/70.0)*e*e*e);
	
	if(p->mpirank==0)
	cout<<"wC 1st: "<<sqrt(9.81*(wH+wdt))<<" | wC Grimshaw: "<<wC;
	

	if(p->mpirank==0)
	cout<<" | wC Fenton: "<<sqrt(9.81*wdt)*sqrt(1.0 + 0.5*e  - (3.0/20.0)*e*e - (3.0/56.0)*e*e*e)<<endl;
	
	X0 = - (2.12*wdt)/(sqrt((0.5*wa)/wdt));
	
	
	if(p->mpirank==0)
	cout<<"X0: "<<X0<<endl;
    
}

void wave_lib_solitary_3rd::wave_prestep(lexer *p, ghostcell *pgc)
{
}
