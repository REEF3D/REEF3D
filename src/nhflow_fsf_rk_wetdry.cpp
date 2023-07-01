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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_rk::wetdry(lexer* p, fdm_nhf* d, ghostcell* pgc, double *UH, double *VH, double *WH, slice &eta)
{
    double wl;
    
    SLICELOOP4
    p->wet_n[IJ] = p->wet[IJ];
    
    
    SLICELOOP4
    {
    //wl =  eta(i,j) + p->wd - d->bed(i,j);
    
        if(d->WL(i,j)>=p->A544)
        p->wet[IJ]=1;
              
        if(d->WL(i,j)<p->A544)
        p->wet[IJ]=0;
      
      
        /*if(p->wet[IJ]==0 && d->eta(i,j)<d->eta(i+1,j))
        p->wet[IJ]=1;
        
        if(p->wet[IJ]==0 && d->eta(i,j)<d->eta(i-1,j))
        p->wet[IJ]=1;
        
        if(p->wet[IJ]==0 && d->eta(i,j)<d->eta(i,j+1))
        p->wet[IJ]=1;
        
        if(p->wet[IJ]==0 && d->eta(i,j)<d->eta(i,j-1))
        p->wet[IJ]=1;*/
        
        //if(d->WL(i,j)<10.0*p->A544)
        //cout<<p->wet[IJ]<<" | "<<wl<<" "<<d->WL(i,j)<<" | "<<p->XP[IP]<<endl;
    }
          

          
    LOOP
    if(p->wet[IJ]==0)
    {
        d->U[IJK] = 0.0;
        d->V[IJK] = 0.0;
        d->W[IJK] = 0.0;
        
        UH[IJK] = 0.0;
        VH[IJK] = 0.0;
        WH[IJK] = 0.0;
    }
    
    FLOOP
    if(p->wet[IJ]==0)
    d->omegaF[FIJK] = 0.0;
    
    
    SLICELOOP4
    if(p->wet[IJ]==0 && p->wet[Im1J]==0 && p->wet[Ip1J]==0 && p->wet[IJm1]==0 && p->wet[IJp1]==0)
    if(eta(i,j)< -p->wd  + d->bed(i,j) - p->A544 + 1.0e-20)
    eta(i,j) = -p->wd  + d->bed(i,j) - p->A544 - 1.0e-15;


    pgc->gcsl_start4(p,d->eta,1);
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    //LOOP
    ///d->test[IJK] = p->wet[IJ];
    
}