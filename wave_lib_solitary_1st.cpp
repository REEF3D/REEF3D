/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wd: "<<wd<<endl;
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
	
	vel = wC * ( ( wH/wd + 3.0*pow(wH/wd,2.0) * (1.0/6.0 - 0.5*pow((wd+z)/wd,2.0)))*(eta/wH)
								-pow(wH/wd,2.0)*(7.0/4.0 - (9.0/4.0)*pow((wd+z)/wd,2.0))*pow(eta/wH,2.0));
								
	if(z<-(p->phimean-p->B126))
	vel=0.0;
											
    return vel;
}

double wave_lib_solitary_1st::wave_w(lexer *p, double x, double y, double z)
{
    double vel,eta;

	eta = wave_eta(p,x,y);
	
	teta = -(wC*(p->simtime) - (x - X0));
	
	vel = wC * sqrt((3.0*wH)/wd)*((wd+z)/wd)*(eta/wd) * tanh(sqrt(0.75*(wH/pow(wd,3.0)))*teta)
		* (1.0 + (wH/(2.0*wd))*(1.0 - 7.0*(eta/wd) - pow((wd+z)/wd,2.0)*(1.0-(3.0*eta)/wH)));
		
	if(z<-(p->phimean-p->B126))
	vel=0.0;
		
    return vel;
}

double wave_lib_solitary_1st::wave_eta(lexer *p, double x, double y)
{
    double eta;
	
	teta = -(wC*(p->simtime) - (x - X0));
	
	eta =  wH/pow(cosh(sqrt(0.75*wH/pow(wd,3.0)) * teta),2.0);
	
	return eta;	
}

double wave_lib_solitary_1st::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_solitary_1st::parameters(lexer *p, ghostcell *pgc)
{
    wH = wa;
	p->wH = wa;
	
	wC = sqrt(9.81*(wH+wd));
	
	X0 = - (2.12*wd)/(sqrt((0.5*wa)/wd));

	
	if(p->mpirank==0)
	cout<<"X0: "<<X0<<endl;
	
	if(p->mpirank==0)
	cout<<"wC: "<<wC<<" wC_old: "<<sqrt(wd*9.81)*(1.0+0.5*(wH/wd))<<endl;
    
}