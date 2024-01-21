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

#include"wave_lib_cnoidal_1st.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_cnoidal_1st::wave_lib_cnoidal_1st(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: 1st-order cnoidal waves "<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<" c: "<<p->wC<<endl;
    cout<<"d/gT^2: "<<wdt/(fabs(p->W22)*wT*wT)<<" H/gT^2: "<<wH/(fabs(p->W22)*wT*wT)<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_cnoidal_1st::~wave_lib_cnoidal_1st()
{
}

double wave_lib_cnoidal_1st::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_cnoidal_1st::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_cnoidal_1st::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel,eta;
	double sn,cn,dn;
	
	teta = 2.0*Km*(p->simtime/wT - x/wL) + pshift;
	
	eta = wave_eta(p,x,y);
	
	elliptic(p,teta,sn,cn,dn);
	
	vel = wC * (eta/wdt - pow(eta/wdt,2.0) + 0.5*(1.0/3.0 - pow((z+wdt)/wdt, 2.0))
	    * wdt*8.0*wH*pow(Km/wL,2.0) * (dn*dn*sn*sn + cn*cn*(-dn*dn + modulus*sn*sn)));

    return vel;
}

double wave_lib_cnoidal_1st::wave_w(lexer *p, double x, double y, double z)
{
    double vel,eta,deta;
	double sn,cn,dn;
	
	teta = 2.0*Km*(p->simtime/wT - x/wL) + pshift;
	
	eta = wave_eta(p,x,y);
	
	elliptic(p,teta,sn,cn,dn);
	
	deta = - wa*sqrt(3.0*wa/wdt)*(1.0/(modulus*wdt))*cn*sqrt(1.0-cn*cn)*sqrt(1.0+modulus*(cn*cn-1.0));
	
	
	vel = - wC * (wdt+z) * (((4.0*wH*(Km/wL)*cn*sn*dn)/wdt) * (1.0 - 2.0*(eta/wdt)) + (1.0/6.0)*wdt
			*(-cn*dn*sn*(dn*dn + modulus*cn*cn - sn*sn) 
			* 64.0*wH*pow(Km/wL,3.0)) * (1.0 - pow((z+wdt)/wdt,2.0)));
			
    return vel;
}

double wave_lib_cnoidal_1st::wave_eta(lexer *p, double x, double y)
{
    double eta;
	double sn,cn,dn;
	
	teta = 2.0*Km*(p->simtime/wT - x/wL);
	
	
	elliptic(p,teta,sn,cn,dn);
	
	eta =  wH*cn*cn + eta2;
	
	return eta;	
}

double wave_lib_cnoidal_1st::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_cnoidal_1st::parameters(lexer *p, ghostcell *pgc)
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

void wave_lib_cnoidal_1st::wave_prestep(lexer *p, ghostcell *pgc)
{
}
