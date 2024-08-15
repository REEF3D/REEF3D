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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"linear_regression_cont.h"
#include<sys/stat.h>
#include<sys/types.h>

void iowave::nhflow_dirichlet_wavegen(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
        double etaval=0.0;
        
        // wave theory
        if(p->B92<20 || p->B92>29)
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        /*d->eta(i-1,j) = (1.0-p->RK_alpha)*eta0(i-1,j)*ramp(p) + p->RK_alpha*eta1(i-1,j)*ramp(p);
        d->eta(i-2,j) = (1.0-p->RK_alpha)*eta0(i-2,j)*ramp(p) + p->RK_alpha*eta1(i-1,j)*ramp(p);
        d->eta(i-3,j) = (1.0-p->RK_alpha)*eta0(i-3,j)*ramp(p) + p->RK_alpha*eta1(i-1,j)*ramp(p);*/
        
        etaval = d->eta(i,j);
        
        d->eta(i-1,j) = etaval;
        d->eta(i-2,j) = etaval;
        d->eta(i-3,j) = etaval;
        
        //cout<<"ETA: "<<eta0(i,j)<<" RK_alpha: "<<p->RK_alpha<<endl;
        }
        
        // wave maker
        if(p->B92>=20 && p->B92<=29)
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        etaval = d->eta(i,j);
        
        d->eta(i-1,j) = etaval;
        d->eta(i-2,j) = etaval;
        d->eta(i-3,j) = etaval;
        }
        
         count=0;
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        

            WETDRYDEEP
            {
            // U, V, W
            uvel = (1.0-p->RK_alpha)*uval0[count] + p->RK_alpha*uval1[count];
            vvel = (1.0-p->RK_alpha)*vval0[count] + p->RK_alpha*vval1[count];
            wvel = (1.0-p->RK_alpha)*wval0[count] + p->RK_alpha*wval1[count];
            
            uvel *= ramp(p);
            vvel *= ramp(p);
            wvel *= ramp(p);
             
            /*uvel = 2.0*uvel - U[IJK];
            vvel = 2.0*vvel - V[IJK];
            wvel = 2.0*wvel - W[IJK];*/
            
                U[Im1JK]=uvel;
                U[Im2JK]=uvel;
                U[Im3JK]=uvel;
                
                V[Im1JK]=vvel;
                V[Im2JK]=vvel;
                V[Im3JK]=vvel;
                
                W[Im1JK]=wvel;
                W[Im2JK]=wvel;
                W[Im3JK]=wvel;
                  
                
            // UH, VH, WH
            uhvel = (1.0-p->RK_alpha)*UHval0[count] + p->RK_alpha*UHval1[count];
            vhvel = (1.0-p->RK_alpha)*VHval0[count] + p->RK_alpha*UHval1[count];
            whvel = (1.0-p->RK_alpha)*WHval0[count] + p->RK_alpha*UHval1[count];
            
            uhvel *= ramp(p);
            vhvel *= ramp(p);
            whvel *= ramp(p);
            
            /*uhvel = 2.0*uhvel - UH[IJK];
            vhvel = 2.0*vhvel - VH[IJK];
            whvel = 2.0*whvel - WH[IJK];*/
            
                UH[Im1JK]=uhvel;
                UH[Im2JK]=uhvel;
                UH[Im3JK]=uhvel;
                
                VH[Im1JK]=vhvel;
                VH[Im2JK]=vhvel;
                VH[Im3JK]=vhvel;
                
                WH[Im1JK]=whvel;
                WH[Im2JK]=whvel;
                WH[Im3JK]=whvel;
            }
            
        ++count;
        }
            
        // ------------------
        if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
		{
            for(int q=0;q<4;++q)
            for(n=0;n<p->gcin_count;++n)
            {
            i=p->gcin[n][0]+q;
            j=p->gcin[n][1];
            k=p->gcin[n][2];
            
            d->EV[IJK]=MIN(d->EV[IJK],1.0e-4);
            }
        pgc->start4V(p,d->EV,24);
		}
        
}
