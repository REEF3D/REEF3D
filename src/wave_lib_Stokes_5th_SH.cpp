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

#include"wave_lib_Stokes_5th_SH.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_Stokes_5th_SH::wave_lib_Stokes_5th_SH(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc)
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: 5th-order Stokes SH waves "<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<endl;
    cout<<"d/gT^2: "<<wdt/(fabs(p->W22)*wT*wT)<<" H/gT^2: "<<wH/(fabs(p->W22)*wT*wT)<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_Stokes_5th_SH::~wave_lib_Stokes_5th_SH()
{
}

double wave_lib_Stokes_5th_SH::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_Stokes_5th_SH::wave_v(lexer *p, double x, double y, double z)
{
    double vel;
	
	vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_Stokes_5th_SH::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = c0*sqrt(fabs(p->W22)/wk)
         *((eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*cos(teta)
         + (pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*cos(2.0*teta)
         + (pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*cos(3.0*teta)
         + (pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*cos(4.0*teta)
         + (pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*cos(5.0*teta));

    return vel;
}

double wave_lib_Stokes_5th_SH::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = c0*sqrt(fabs(p->W22)/wk)
         *((eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*sinh(wk*(wdt+z))*sin(teta)
         + (pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*sinh(2.0*wk*(wdt+z))*sin(2.0*teta)
         + (pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*sinh(3.0*wk*(wdt+z))*sin(3.0*teta)
         + (pow(eps,4.0)*a44)*sinh(4.0*wk*(wdt+z))*sin(4.0*teta)
         + (pow(eps,5.0)*a55)*sinh(5.0*wk*(wdt+z))*sin(5.0*teta));

    return vel;
}

double wave_lib_Stokes_5th_SH::wave_eta(lexer *p, double x, double y)
{
    double eta;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    eta =  (1.0/wk)*((eps + pow(eps,3.0)*b31 - pow(eps,5.0)*(b53 + b55))*cos(teta)
                    + (pow(eps,2.0)*b22 + pow(eps,4.0)*b42)*cos(2.0*teta)
                    + (-pow(eps,3.0)*b31 + pow(eps,5.0)*b53)*cos(3.0*teta)
                    + pow(eps,4.0)*b44*cos(4.0*teta)
                    + pow(eps,5.0)*b55*cos(5.0*teta));

    return eta;
}

double wave_lib_Stokes_5th_SH::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_Stokes_5th_SH::parameters(lexer *p, ghostcell *pgc)
{
    eps = 0.5*wk*wH;


    S = 1.0/cosh(2*wk*wdt);
    C = 1.0 - S;


    a11 = 1.0/sinh(wk*wdt);

    a22 = 3.0*S*S/(2.0*C*C);

    a31 = (-4.0 - 20.0*S + 10.0*S*S -13.0*S*S*S)/(8.0*sinh(wk*wdt)*C*C*C);

    a33 = (-2.0*S*S + 11*S*S*S)/(8.0*sinh(wk*wdt)*C*C*C);

    a42 = (12.0*S - 14.0*S*S - 264.0*S*S*S - 45.0*pow(S,4.0) - 13.0*pow(S,5.0))/(24.0*pow(C,5.0));

    a44 = (10.0*S*S*S - 174.0*pow(S,4.0) + 291.0*pow(S,5.0) + 278.0*pow(S,6.0))/(48.0*(3.0 + 2.0*S)*pow(C,5.0));

    a51 = (-1184.0 + 32.0*S + 13232.0*S*S + 21712.0*S*S*S + 20940.0*pow(S,4.0) + 12554.0*pow(S,5.0)
           -500.0*pow(S,6.0) - 3341.0*pow(S,7.0) - 670.0*pow(S,8.0))/(64.0*sinh(wk*wdt)*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));

    a53 = (4.0*S + 105.0*pow(S,2.0) + 198.0*pow(S,3.0) - 1376.0*pow(S,4.0) - 1302.0*pow(S,5.0) - 117.0*pow(S,6.0) + 58.0*pow(S,7.0))
          /(32.0*sinh(wk*wdt)*(3.0 + 2.0*S)*pow(C,6.0));

    a55 = (-6.0*S*S*S + 272.0*pow(S,4.0) - 1552.0*pow(S,5.0) + 852.0*pow(S,6.0) + 2029.0*pow(S,7.0) + 430.0*pow(S,8.0))
          /(64.0*sinh(wk*wdt)*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));


    b22 = ((cosh(wk*wdt)/sinh(wk*wdt))*(1.0 + 2.0*S))/(2.0*C);

    b31 = (-3.0*(1.0 + 3.0*S + 3.0*S*S + 2.0*S*S*S))/(8.0*C*C*C);

    b42 = ((cosh(wk*wdt)/sinh(wk*wdt))*(6.0 - 26.0*S - 182.0*S*S - 204.0*S*S*S - 25.0*pow(S,4.0) + 26*pow(S,5.0)))
          /(6.0*(3.0 + 2.0*S)*pow(C,4.0));

    b44 = ((cosh(wk*wdt)/sinh(wk*wdt))*(24.0 + 92.0*S + 122.0*S*S + 66.0*S*S*S + 67.0*pow(S,4.0) + 34.0*pow(S,5.0)))
          /(24.0*(3.0 + 2.0*S)*pow(C,4.0));

    b53 = (9.0*(132.0 + 17.0*S - 2216.0*S*S - 5897.0*S*S*S - 6292.0*pow(S,4.0) - 2687.0*pow(S,5.0) + 194.0*pow(S,6.0)
                + 467.0*pow(S,7.0) + 82.0*pow(S,8.0)))/(128.0*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));

    b55 = (5.0*(300.0 + 1579.0*S + 3176.0*S*S + 2949.0*S*S*S + 1188.0*pow(S,4.0) + 675.0*pow(S,5.0) + 1326.0*pow(S,6.0)
                + 827.0*pow(S,7.0) + 130.0*pow(S,8.0)))/(384.0*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));


    c0 = sqrt(tanh(wk*wdt));

    c2 = (c0*(2.0 + 7.0*S*S)/(4.0*C*C));

    c4 = (c0*(4.0 + 32.0*S -116.0*S*S - 400.0*S*S*S - 71.0*pow(S,4.0) + 146.0*pow(S,5.0)))/(32.0*pow(C,5.0));


    e2 = (tanh(wk*wdt)*(2.0 + 2.0*S + 5.0*S*S))/(4.0*pow(1.0 - S,2.0));

    e4 = (tanh(wk*wdt)*(8.0 + 12.0*S - 152.0*S*S - 308.0*pow(S,3.0) - 42.0*pow(S,4.0) + 77.0*pow(S,5.0)))/(32.0*pow(1.0 - S, 5.0));
}

void wave_lib_Stokes_5th_SH::wave_prestep(lexer *p, ghostcell *pgc)
{
}
