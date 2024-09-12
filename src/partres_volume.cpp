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

double partres::volume(lexer *p, fdm &a, particles_obj &PP)
{
        double sum=0;
        ILOOP
            JLOOP
                sum += columnSum[IJ];

        return PI*pow(PP.d50,3.0)*(sum)/6.0;
}

    /// @brief Calculate number of particles in cell ( \p i , \p j , \p k )
void partres::particlePerCell(lexer *p, ghostcell &pgc, particles_obj &PP)
{
        PLAINLOOP
        cellSum[IJK]=0;

        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]>INT32_MIN)
            {
                i=p->posc_i(PP.X[n]);
                j=p->posc_j(PP.Y[n]);
                k=p->posc_k(PP.Z[n]);
                cellSum[IJK] += PP.ParcelFactor[n];
            }
        
        pgc.start4V_par(p,cellSum,11);
        pgc.start4V_par(p,cellSumTopo,11);
}

    /// @brief Topo volume in cell div. by particle volume
    /// Uses i,j&k from increment to pass cell identifier
    /// @param d50 Sauter diameter of particles
    /// @return Ceil of number of particles in cell IJK
double partres::maxParticlesPerCell(lexer *p, fdm &a, double d50, bool topo, bool cell)
{   
        double DZN=topo?0:p->DZN[KP];

        if(topo)
        {
            if (a.topo(i,j,k)<=-0.5*p->DZN[KP]+1.0e-13)
            DZN=p->DZN[KP];
            else if(a.topo(i,j,k)<0.5*p->DZN[KP] -5.0e-18)
            DZN=(p->DZN[KP]*0.5 + a.topo(i,j,k));
        }

        return 6.0*p->DXN[IP]*p->DYN[JP]*DZN*(1.0+(cell?0:-a.porosity(i,j,k)))/(PI*pow(d50,3.0));
}
