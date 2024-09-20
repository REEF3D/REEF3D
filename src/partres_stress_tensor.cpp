/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

    /// @brief Calculate intra-particle stress trensor for cell ( \p i , \p j , \p k )
void partres::ParticleStressTensor(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
{
    BLOOP
    {
    Ps = 5.0;
    beta = 3.5;
    epsilon = 10e-7;
    Tc = 0.6;
    
    Ts = PI*pow(PP.d50,3.0)*(cellSum[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
    
    Ts = MAX(Ts,0.0);
    Ts = MIN(Ts,1.0);

    stressTensor[IJK] = Ps*pow(Ts,beta)/MAX(Tc-Ts,epsilon*(1.0-Ts));
    
    a.test(i,j,k) = stressTensor[IJK];
    
    //cout<<stressTensor[IJK]<<" "<<Ts<<endl;
    }
    
    pgc.start4V_par(p,stressTensor,10);
}

    /// @brief Calculate solid volume fraction for cell ( \p i , \p j , \p k )
double partres::theta_s(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k) const
{   
        double theta = PI*pow(PP.d50,3.0)*(cellSum[IJK]+cellSumTopo[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
        if(theta>1)
        theta=1;
        if(theta<0)
        theta=0;
        return theta;
}    

