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
        
        if(p->A515==1)
        etaval = d->eta(i,j);
        
        if(p->A515==2)
        etaval = eta(i,j);
        
        d->eta(i-1,j) = etaval;
        d->eta(i-2,j) = etaval;
        d->eta(i-3,j) = etaval;
        }
        
        // wave maker
        if(p->B92>=20 && p->B92<=29)
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        if(p->A515==1)
        etaval = d->eta(i,j);
        
        if(p->A515==2)
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
        

            //WETDRYDEEP
            //{
            // U, V, W
            uvel = uval[count];
            vvel = vval[count];
            wvel = wval[count];
            
            uvel *= ramp(p);
            vvel *= ramp(p);
            wvel *= ramp(p);
             
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
            uhvel = UHval[count];
            vhvel = VHval[count];
            whvel = WHval[count];
            
            uhvel *= ramp(p);
            vhvel *= ramp(p);
            whvel *= ramp(p);
            
                UH[Im1JK]=uhvel;
                UH[Im2JK]=uhvel;
                UH[Im3JK]=uhvel;
                
                VH[Im1JK]=vhvel;
                VH[Im2JK]=vhvel;
                VH[Im3JK]=vhvel;
                
                WH[Im1JK]=whvel;
                WH[Im2JK]=whvel;
                WH[Im3JK]=whvel;
            //}
            
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
        pgc->start24V(p,d->EV,24);
        }
        
}
