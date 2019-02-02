/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"wave_lib_Stokes_2nd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_Stokes_2nd::wave_lib_Stokes_2nd(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc)
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: 2nd-order Stokes waves; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wd: "<<wd<<" kd: "<<wd*wk<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_Stokes_2nd::~wave_lib_Stokes_2nd()
{
}

double wave_lib_Stokes_2nd::wave_u(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);
    
    return cosgamma*vel;
}

double wave_lib_Stokes_2nd::wave_v(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_Stokes_2nd::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = ww*wa*( sinh(wk*(wd+z))/sinh(wk*(wd)) ) * sin(teta)
         + 0.75*wk*ww*wa*wa*( sinh(2.0*wk*(wd+z))/pow(sinh(wk*(wd)),4.0) ) * sin(2.0*teta);

    return vel;
}

double wave_lib_Stokes_2nd::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = ww*wa*( cosh(wk*(wd+z))/sinh(wk*(wd)) ) * cos(teta)
         + 0.75*wk*ww*wa*wa*( cosh(2.0*wk*(wd+z))/pow(sinh(wk*(wd)),4.0) ) * cos(2.0*teta);

    return vel;
}

double wave_lib_Stokes_2nd::wave_eta(lexer *p, double x, double y)
{
    double eta;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    eta =  wa*cos(teta) + 0.25*wk*wa*wa*(cosh(wk*wd)/pow(sinh(wk*wd),3.0)) * (2.0 + cosh(2.0*wk*wd)) * cos(2.0*teta);

    return eta;
}

double wave_lib_Stokes_2nd::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
	teta = wk*x-ww*(p->simtime) + pshift;

    fi = ((ww*wa)/wk)*( cosh(wk*(wd+z))/sinh(wk*(wd)) ) * sin(teta)
    
         + (3.0/8.0)*ww*wa*wa*( cosh(2.0*wk*(wd+z))/pow(sinh(wk*(wd)),4.0) ) * sin(2.0*teta);


    return fi;
}

void wave_lib_Stokes_2nd::parameters(lexer *p, ghostcell *pgc)
{
    
}