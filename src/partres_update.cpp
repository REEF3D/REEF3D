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

    /**
     * @brief Updates the topo field.
     *
     * This function is responsible for updating the topo field.
     * It does so by calculating the difference between the columnSum and the count value.
     *
     * @param p A pointer to the lexer object.
     * @param pgc A reference to the ghostcell object.
     * @param topo A reference to the field4a object.
     * @param d50 The diameter value.
     */
void partres::update(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
{
    int count=0;
    ILOOP
    JLOOP
    {
    if(bedChange[IJ]!=0)
    {
        KLOOP
        {
            a.topo(i,j,k) -= bedChange[IJ]*1.0/6.0*PI*pow(PP.d50,3)/(p->DXN[IP]*p->DYN[JP]);
            a.fb(i,j,k) = bedChange[IJ];
        }
        columnSum[IJ] += bedChange[IJ];
        if(bedChange[IJ]<0)
        count+=activateNew(p,a,PP);
        bedChange[IJ] = 0;
    }
    }
    pgc.start4a(p,a.topo,150);
    if(count>0)
    cout<<"On partion "<<p->mpirank<<" were "<<count<<" additional particles activated."<<endl;

        // double count;
        // ILOOP
        // {
        //     if(p->XN[IP]>=p->global_xmin+p->Q73)
        //         JLOOP
        //         {
        //             count = 0.0;
        //             KLOOP
        //             {
        //                 count += cellSum[IJK] + cellSumTopo[IJK];
        //                 if(k>0 && cellSumTopo[IJKm1]==0)
        //                 break;
        //             }
        //             if(count != columnSum[IJ])
        //             {
        //                 KLOOP
        //                 topo(i,j,k) -= (count-columnSum[IJ])*4.0/3.0*PI*pow(d50/2.0,3)/(p->DXN[IP]*p->DYN[JP]);
        //             }
        //             columnSum[IJ] = count;
        //         }
        // }
        // pgc.start4a(p,topo,150);
}