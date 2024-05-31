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

void iowave::nhflow_dirichlet_wavegen(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    if(p->count<=1)
    {
    netQ = 0.0;
    netQn = 0.0;
    netQnn = 0.0;
    netQnnn = 0.0;
    netQ_corr = 0.0;
    netQ_minpeak=0.0;
    netQ_maxpeak=0.0;
    netQ_minpeak_n=0.0;
    netQ_maxpeak_n=0.0;
    netQ_minpeak_nn=0.0;
    netQ_maxpeak_nn=0.0;
    }
    
    netQnnn = netQnn;
    netQnn = netQn;
    netQn = netQ;
    
    netQ_corr = -0.001*netQ_avgpeak/p->dt;
    
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        d->eta(i-1,j) = d->eta(i,j);
        d->eta(i-2,j) = d->eta(i,j);
        d->eta(i-3,j) = d->eta(i,j);
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
            uvel=uval[count];// + netQ_corr/d->WL(i,j);
            vvel=vval[count];
            wvel=wval[count];
            
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
            uvel=UHval[count] + netQ_corr;
            vvel=VHval[count];
            wvel=WHval[count];
            
            uvel *= ramp(p);
            vvel *= ramp(p);
            wvel *= ramp(p);
            
            /*uvel = 2.0*uvel - UH[IJK];
            vvel = 2.0*vvel - VH[IJK];
            wvel = 2.0*wvel - WH[IJK];*/
            
                UH[Im1JK]=uvel;
                UH[Im2JK]=uvel;
                UH[Im3JK]=uvel;
                
                VH[Im1JK]=vvel;
                VH[Im2JK]=vvel;
                VH[Im3JK]=vvel;
                
                WH[Im1JK]=wvel;
                WH[Im2JK]=wvel;
                WH[Im3JK]=wvel;
                
                netQ += p->dt*uvel;
            }
            
        ++count;
		}
        
        // ------------------
        // netQ
        netQ_max = MAX(netQ_max,netQ);
        netQ_min = MIN(netQ_min,netQ);
        
        if(netQnn>netQnnn && netQn>netQnn && netQ<netQn)
        {
        netQ_maxpeak_nn = netQ_maxpeak_n;
        netQ_maxpeak_n = netQ_maxpeak;
        netQ_maxpeak = netQn;
        }
        
        if(netQnn<netQnnn && netQn<netQnn && netQ>netQn)
        {
        netQ_minpeak_nn = netQ_minpeak_n;
        netQ_minpeak_n = netQ_minpeak;
        netQ_minpeak = netQn;
        }
        
        
        netQ_avgpeak = (1.0/6.0)*(1.0*netQ_minpeak + 1.0*netQ_maxpeak + 1.0*netQ_minpeak_n + 1.0*netQ_maxpeak_n + 1.0*netQ_minpeak_nn + 1.0*netQ_maxpeak_nn);
        
        //netQ_avgpeak = 0.5*(netQ_minpeak + netQ_maxpeak);
            
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
        
        //netQ=pgc->globalsum(netQ);
        
        if(p->mpirank==0)
        {
        cout<<"netQ: "<<netQ<<" netQ_maxpeak: "<<netQ_maxpeak<<" netQ_minpeak: "<<netQ_minpeak<<endl;
        cout<<"netQ_max: "<<netQ_max<<" netQ_min: "<<netQ_min<<endl;
        cout<<"netQ_avgpeak: "<<netQ_avgpeak<<" netQ_corr: "<<netQ_corr<<endl;
        }
        
        
}
