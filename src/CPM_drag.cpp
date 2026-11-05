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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

// drag per unit particle mass Dp [1/s], force = m Dp (Uf - Up); vel: |Uf - Up|
//
// Q 51 1: Andrews & O'Rourke (1996)
//   Dp = Cd 3/8 rho_f/rho_p |Uf-Up|/r_p,  Cd = 24/Re (thf^-2.65 + 1/6 Re^2/3 thf^-1.78)
//   written without the 1/Re singularity:  Dp = 18 nu rho_f/(rho_p d^2) (thf^-2.65 + 1/6 Re^2/3 thf^-1.78)
//
// Q 51 2: Gidaspow (1994), the standard of two-fluid and MP-PIC models for dense beds
//   thf < 0.8: Ergun        Dp = 150 ths mu/(thf rho_p d^2) + 1.75 rho_f |Uf-Up|/(rho_p d)
//   thf >= 0.8: Wen & Yu    Dp = 18 mu/(rho_p d^2) (1 + 0.15 Re^0.687) thf^-2.65,  Re = thf rho_f |Uf-Up| d/mu
double CPM::drag_model(lexer *p, double d, double rhoS, double vel, double Ts)
{
    double Tf = MAX(1.0-Ts, 1.0-theta_max);
    Tf = MIN(Tf,1.0);
    
    vel = fabs(vel);
    
    double mu = p->W2*p->W1;
    double Dp;
    
    if(p->Q51==2)
    {
        if(Tf<0.8)
        Dp = 150.0*(1.0-Tf)*mu/(Tf*rhoS*d*d) + 1.75*p->W1*vel/(rhoS*d);
        
        else
        {
            double Re = Tf*p->W1*vel*d/mu;
            Dp = 18.0*mu/(rhoS*d*d) * (Re<1000.0 ? 1.0 + 0.15*pow(Re,0.687) : 0.44*Re/24.0) * pow(Tf,-2.65);
        }
    }
    
    else
    {
        double Rep = vel*d/p->W2;
        Dp = 18.0*p->W2*p->W1/(rhoS*d*d) * (pow(Tf,-2.65) + (1.0/6.0)*pow(Rep,2.0/3.0)*pow(Tf,-1.78));
    }

    return Dp;
}
