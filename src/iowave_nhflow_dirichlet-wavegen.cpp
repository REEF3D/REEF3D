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
    
    netQ = 0.0;
    
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        d->eta(i-1,j) = eta(i,j)*ramp(p);
        d->eta(i-2,j) = eta(i,j)*ramp(p);
        d->eta(i-3,j) = eta(i,j)*ramp(p);
        }
        
         count=0;
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        
            WETDRYDEEP
            {
            double uhcorr = 0.0;
            double ucorr = 0.0;
            
            /*
            if(p->count>2);
            {
            uhcorr = netV_corr/(p->DYN[JP]*p->dt);
            ucorr /= (d->eta(i,j)+d->depth(i,j));
            }*/
            
            uhcorr = 0.0;
            ucorr = 0.0;
            
            // U, V, W
            uvel=uval[count];
            vvel=vval[count];
            wvel=wval[count];
            
            uvel *= ramp(p);
            vvel *= ramp(p);
            wvel *= ramp(p);
            
                U[Im1JK]=uvel + ucorr;
                U[Im2JK]=uvel + ucorr;
                U[Im3JK]=uvel + ucorr;
                
                V[Im1JK]=vvel;
                V[Im2JK]=vvel;
                V[Im3JK]=vvel;
                
                W[Im1JK]=wvel;
                W[Im2JK]=wvel;
                W[Im3JK]=wvel;
                
                
            // UH, VH, WH
            uvel=UHval[count];
            vvel=VHval[count];
            wvel=WHval[count];
            
            uvel *= ramp(p);
            vvel *= ramp(p);
            wvel *= ramp(p);
            
                UH[Im1JK]=uvel + uhcorr;
                UH[Im2JK]=uvel + uhcorr;
                UH[Im3JK]=uvel + uhcorr;
                
                VH[Im1JK]=vvel;
                VH[Im2JK]=vvel;
                VH[Im3JK]=vvel;
                
                WH[Im1JK]=wvel;
                WH[Im2JK]=wvel;
                WH[Im3JK]=wvel;
                
                netQ += (uvel*p->DYN[JP] + uhcorr);
            }
            
        ++count;
		}
        
        // ------------------
        // netQ
        double relaxfac = 0.0005;
        
        netQ = pgc->globalsum(netQ);
        
        netV += p->dt*0.5*(3.0*netQ - netQ_n)*ramp_corr(p);
        
        linreg->linreg_cont_func(p,pgc,p->simtime,netV,b0,b1);

        
        if(p->count>1)
        {
        netV_corr_n = netV_corr;
        netV_corr =  -relaxfac*(b1*(p->simtime) + b0)*ramp_corr(p);
        /*
        if(netV_corr>netV_corr_n)
        netV_corr = MAX(netV_corr, netV_corr + 0.1*(netV_corr-netV_corr_n));
        
        if(netV_corr<netV_corr_n)
        netV_corr = MIN(netV_corr, netV_corr + 0.1*(netV_corr-netV_corr_n));*/
        }
        
        
        netQ_n = netQ;
        
        /* 
        UH = WL*U
        
        A = WL*DY
        
        Q = A*U
        V = Q*dt
        V = A*U*dt
        
        U_corr = dV/(A*dt)
        UH_corr = U*WL
        UH_corr = dV*WL/(WL*DY*dt)
        UH_corr = dV/(DY*dt)
        */
        
        /*
        
        or: correct via eta....
        
        */
            
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
        
        
        
        if(p->mpirank==0)
        {
        cout<<"netV: "<<netV<<" netQ: "<<netQ<<" netV_corr: "<<netV_corr<<endl;
        cout<<setprecision(6)<<"netQ_b0: "<<setprecision(6)<<b0<<" netQ_b1: "<<b1<<endl;
        }
        
    if(p->mpirank==0 && p->count%p->P12==0)
	 {
     logout<<p->simtime+p->dt<<"\t"<<netV<<endl; 
	 }
    
}
