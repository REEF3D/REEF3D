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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

// Snider (2001): Ps theta^beta / max(theta_cs - theta, eps (1-theta))
void CPM::stress_snider(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const double Ps = p->Q14;
    const double beta = p->Q15;
    const double eps = p->Q16;
    const double thcs = p->Q17;
    
    double denom,dPdT;
    
    cmax=0.0;

    BASELOOP
    {
        denom = MAX(thcs-Ts(i,j,k), eps*(1.0-Ts(i,j,k)));
        
        Tau(i,j,k) = Ps*pow(MAX(Ts(i,j,k),0.0),beta)/denom;
        
        // wave speed of the particle stress, dTau/dtheta
        dPdT = Tau(i,j,k)*(beta/MAX(Ts(i,j,k),1.0e-6) + (thcs-Ts(i,j,k)>eps*(1.0-Ts(i,j,k))?1.0/denom:0.0));
        
        cmax = MAX(cmax,dPdT);
    }
    
    cmax = sqrt(pgc->globalmax(cmax)/p->S22);

    pgc->start4a(p,Tau,1);
}
