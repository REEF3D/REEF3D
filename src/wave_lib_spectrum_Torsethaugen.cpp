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

#include"wave_lib_spectrum.h"
#include"lexer.h"
#include"ghostcell.h"

double wave_lib_spectrum::Torsethaugen(lexer *p, double w)
{
	double Tpf,Tp1,Tp2,H1,H2,Rs,Rw,f,f1,f2;
	double af,ae,au,a10,a1,kg,b1,a20,a2,a3;
	double gamma,gammaf,Agamma,sigma,sp,epsl,epsu,Tl,Tu,G0;
	double E1,E2,S1,S2;
	
	af = 6.6;
	ae = 2.0;
	au = 25.0;
	a10 = 0.7;
	a1 = 0.5;
	kg  = 35.0;
	b1 = 2.0;
	a20 = 0.6;
	a2 = 0.3;
	a3 = 6.0;
	G0 = 3.26;
	
	f = w/(2.0*PI);
	
	Tpf = af*pow(p->wHs,1.0/3.0);
	Tl = ae*sqrt(p->wHs);
	Tu = au;
	
	epsl = (Tpf-p->wTp)/(Tpf-Tl);
	epsu = (p->wTp-Tpf)/(Tu-Tpf);
	
	
// Wind dominated sea
	if(p->wTp<=Tpf)
	{
	// Primary Peak
	Rw = (1.0-a10)*exp(-pow(epsl/a1,2.0)) + a10;
	H1 = Rw*p->wHs;	
	
	Tp1 = p->wTp;	
	
	sp = ((2.0*PI)/9.81)*(H1/pow(Tp1,2.0));
	gamma = kg*pow(sp,(6.0/7.0));
	
	// Secondary Peak
	H2 = sqrt(1.0 - Rw*Rw)*p->wHs;
	
	Tp2 = Tpf + b1;
	}
	
// Swell dominaded sea
	if(p->wTp>Tpf)
	{
	// Primary Peak
	Rs = (1.0-a20)*exp(-pow(epsu/a2,2.0)) + a20;
	H1 = Rs*p->wHs;	
	
	Tp1 = p->wTp;	
	
	sp = ((2.0*PI)/9.81)*(p->wHs/pow(Tpf,2.0));
	gammaf = kg*pow(sp,(6.0/7.0));
	gamma = gammaf*(1.0 + a3*epsu);
	
	// Secondary Peak
	H2 = sqrt(1.0 - Rs*Rs)*p->wHs;
	
	Tp2 = af*pow(H2,(1.0/3.0));
	}
	
    E1 = (1.0/16.0)*pow(H1,2.0)*Tp1;
	E2 = (1.0/16.0)*pow(H2,2.0)*Tp2;
	
	gamma=MAX(1.0,gamma);
	Agamma = (1.0 + 1.1*pow(log(gamma),1.19))/gamma;
	
	
	
	f1 = f*Tp1;
	f2 = f*Tp2;
	
	f1 = fabs(f1)>1.0e-20?f1:1.0e20;
	f2 = fabs(f2)>1.0e-20?f2:1.0e20;
	
	if(f1<=1.0)
	sigma=0.07;
	
	if(f1>1.0)
	sigma=0.09;
	
	S1 = G0*Agamma * pow(f1,-4.0) * exp(-pow(f1,-4.0)) * pow(gamma,exp(- (1.0/(2.0*pow(sigma,2.0)))*pow(f1-1.0,2.0) ));
	
	S2 = G0*pow(f2,-4.0) * exp(-pow(f2,-4.0));
	
	//cout<<Tp1<<" "<<Tp2<<" | "<<S1<<" "<<S2<<" | "<<Tpf<<endl;
	
	Sval = (E1*S1 + E2*S2 )/(2.0*PI);
	
    return Sval;
}




