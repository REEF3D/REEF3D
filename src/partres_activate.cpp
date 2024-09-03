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

#include "partres.h"
#include "particles_obj.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

int partres::activateNew(lexer *p, fdm &a, particles_obj &PP)
{
        double tolerance = 5e-18;
        double x,y,z,ipolTopo,ipolSolid;
        int flag=0;
        size_t index;

        double counter=0;
        int maxTries=1000;
        int tries=0;
        int count=0;

        if(PP.size-bedChange[IJ]>0.9*PP.capacity)
            PP.reserve();

        while(counter<-bedChange[IJ] && tries<maxTries)
        {   
            x = p->XN[IP] + p->DXN[IP]*double(rand() % 10000)/10000.0;
            y = p->YN[JP] + p->DYN[JP]*double(rand() % 10000)/10000.0;
            k = 0;
            z = p->ZN[KP]+0.5*p->DZP[KP]-a.topo(i,j,k) - 5.0*PP.d50*double(rand() % 10000)/10000.0;
            k = p->posc_k(z);

            ipolTopo = p->ccipol4_b(a.topo,x,y,z);
            ipolSolid = p->ccipol4_b(a.solid,x,y,z);

            if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0))
                if(cellSumTopo[IJK]>=p->Q41)
                {
                    index = PP.add(x,y,z,flag,0,0,0,p->Q41);
                    counter += PP.PackingFactor[index];
                    cellSumTopo[IJK] -= PP.PackingFactor[index];
                    ++count;
                }
                else if (cellSumTopo[IJK]+cellSumTopo[IJKm1]>=p->Q41)
                {
                    index = PP.add(x,y,z,flag,0,0,0,p->Q41);
                    counter += PP.PackingFactor[index];
                    cellSumTopo[IJK] -= PP.PackingFactor[index];
                    cellSumTopo[IJKm1] += cellSumTopo[IJK];
                    cellSumTopo[IJK] = 0;

                    ++count;
                    break;
                }
            
            ++tries;
        }
        return count;
}
