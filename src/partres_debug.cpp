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
#include"field4a.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void partres::debug(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s)
{
    /*
        PLAINLOOP
        {
            a.test(i,j,k) = cellSumTopo[IJK];
            a.vof(i,j,k) = (s.tau_eff[IJ]>s.tau_crit[IJ]?-1:1);
        }
        double topoDist=0;
        double u=0,v=0;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]>INT32_MIN)
            {
                if(PP.Flag[n]==0)
                {
                    // PP.shear_eff[n]=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                    // PP.shear_crit[n]=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                }

                topoDist=p->ccipol4_a(a.topo,PP.X[n],PP.Y[n],PP.Z[n]);

                if(topoDist<velDist*p->DZP[KP])
                {
                    u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                    v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                }
                else
                {
                    u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]);
                    v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]);
                }

                PP.Uf[n]=u;
                PP.Vf[n]=v;
            }*/
}