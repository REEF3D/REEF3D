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

#include"wave_lib_solitary_1st.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_solitary_1st::wave_lib_solitary_1st(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: 1st-order solitary wave; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_solitary_1st::~wave_lib_solitary_1st()
{
}

double wave_lib_solitary_1st::wave_u(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);
											
    return cosgamma*vel;
}

double wave_lib_solitary_1st::wave_v(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);
											
    return singamma*vel;
}

double wave_lib_solitary_1st::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel,eta;

	eta = wave_eta(p,x,y);
	
	vel = wC * ( ( wH/wdt + 3.0*pow(wH/wdt,2.0) * (1.0/6.0 - 0.5*pow((wdt+z)/wdt,2.0)))*(eta/wH)
								-pow(wH/wdt,2.0)*(7.0/4.0 - (9.0/4.0)*pow((wdt+z)/wdt,2.0))*pow(eta/wH,2.0));
                            
    return vel;
}

double wave_lib_solitary_1st::wave_w(lexer *p, double x, double y, double z)
{
    double vel,eta;

	eta = wave_eta(p,x,y);
	
	teta = -(wC*(p->simtime) - (x - X0));
	
	vel = wC * sqrt((3.0*wH)/wdt)*((wdt+z)/wdt)*(eta/wdt) * tanh(sqrt(0.75*(wH/pow(wdt,3.0)))*teta)
		* (1.0 + (wH/(2.0*wdt))*(1.0 - 7.0*(eta/wdt) - pow((wdt+z)/wdt,2.0)*(1.0-(3.0*eta)/wH)));
		
    return vel;
}

double wave_lib_solitary_1st::wave_eta(lexer *p, double x, double y)
{
    double eta;
	
	teta = -(wC*(p->simtime) - (x - X0));
	
	eta =  wH/pow(cosh(sqrt(0.75*wH/pow(wdt,3.0)) * teta),2.0);
	
	return eta;	
}

double wave_lib_solitary_1st::wave_fi(lexer *p, double x, double y, double z)
{
    double fi=0.0;
    
    return fi;
}

void wave_lib_solitary_1st::parameters(lexer *p, ghostcell *pgc)
{	
	wC = sqrt(9.81*(wH+wdt));
	
	X0 = - (2.12*wdt)/(sqrt((0.5*wa)/wdt));

	
	if(p->mpirank==0)
	cout<<"X0: "<<X0<<endl;
	
	if(p->mpirank==0)
	cout<<"wC: "<<wC<<" wC_old: "<<sqrt(wdt*9.81)*(1.0+0.5*(wH/wdt))<<endl;
    
}

void wave_lib_solitary_1st::wave_prestep(lexer *p, ghostcell *pgc)
{
}
