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

#include"wave_lib_cnoidal_shallow.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_cnoidal_shallow::wave_lib_cnoidal_shallow(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: shallow cnoidal waves; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_cnoidal_shallow::~wave_lib_cnoidal_shallow()
{
}

double wave_lib_cnoidal_shallow::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_cnoidal_shallow::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_cnoidal_shallow::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel,eta;
	
	eta = wave_eta(p,x,y);

    vel = sqrt(9.81/wdt) * eta; 

    return vel;
}

double wave_lib_cnoidal_shallow::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
	double sn,cn,dn;
	
	teta = 2.0*Km*(x/wL - p->simtime/wT) + pshift;

	elliptic(p,teta,sn,cn,dn);

    vel = sqrt(9.81/wdt)*wH * ((4.0*Km*wdt)/wL) * ((wdt+z)/wdt) * cn*sn*dn;

    return vel;
}

double wave_lib_cnoidal_shallow::wave_eta(lexer *p, double x, double y)
{
    double eta;
	double sn,cn,dn;
	
	teta = 2.0*Km*(x/wL - p->simtime/wT) + pshift;
	
	elliptic(p,teta,sn,cn,dn);
	
	eta =  wH*cn*cn + eta2;
	
	return eta;	
}

double wave_lib_cnoidal_shallow::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_cnoidal_shallow::parameters(lexer *p, ghostcell *pgc)
{
    double diff=1.0;
	int qq,maxiter;
	double modulus_old=0.5;
	double Ur;
	modulus = 0.9;
	maxiter =5000;
	
	Ur = (wH*wL*wL)/pow(wdt,3.0);
	if(p->mpirank==0)
	cout<<"Ursell number: "<<Ur<<endl;
	
	qq=1;
	while(fabs(diff)>1.0e-12)
	{
		++qq;
		modulus_old = modulus;
		
		
		modulus = sqrt((3.0/16.0)* ((wH*wL*wL)/(wdt*wdt*wdt*modulus_old*pow(K_elliptic_1(modulus_old),2.0))));
		
		diff = modulus - modulus_old;
		
		modulus = modulus_old + 0.1*diff;
		
		modulus=MIN(0.999, modulus);
	
		if(qq>maxiter)
		break;
	}
	//modulus=MAX(0.9, modulus);
	if(p->mpirank==0)
	cout<<"MODULUS: "<<modulus<<"    qq: "<<qq<<endl;	
	
	Km = K_elliptic_1(modulus);
	Em = E_elliptic_1(modulus);
	
	
	eta2 = (1.0/modulus - 1.0 - Em/(Km*modulus))*wH;
	
	//wC = sqrt((9.81*wdt)/(1.0 + wH/wdt*(1.0/pow(modulus,2.0) -2.0)));
	
	wC = sqrt(9.81*wdt*(1.0 +(wH/wdt)*(2.0/modulus - 1.0 - 3.0/modulus*Em/Km)));
	
	if(p->mpirank==0)	
	{
	cout<<"WAVE TROUGH: "<<eta2+wdt<<endl;
	cout<<"wC: "<<wC<<" wC_old: "<<(wL/wT)<<endl;
	}
    
}

void wave_lib_cnoidal_shallow::wave_prestep(lexer *p, ghostcell *pgc)
{
}
