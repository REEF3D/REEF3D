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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_active_wavegen(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W)
{
    double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
    
        // wavegen
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];		

        uvel=uval[count]*ramp(p);
        vvel=vval[count]*ramp(p);
        wvel=wval[count]*ramp(p);

        U[Im1JK]=uvel+p->Ui;
        U[Im2JK]=uvel+p->Ui;
        U[Im3JK]=uvel+p->Ui;
            
        V[Im1JK]=vvel;
        V[Im2JK]=vvel;
        V[Im3JK]=vvel;
            
        W[Im1JK]=wvel;
        W[Im2JK]=wvel;
        W[Im3JK]=wvel;
            
                
                // fsf deviation
                eta_T = wave_eta(p,pgc,x,0.0);
                eta_M = d->eta(i,j); 
                eta_R = eta_T-eta_M;
				
                if(p->B98==4)
                Uc=eta_R*sqrt(9.81/p->wd);
                
                U[Im1JK] += Uc;
                U[Im2JK] += Uc;
                U[Im3JK] += Uc;
 
        ++count;
		}
        
        
         if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
		{
            for(int q=0;q<4;++q)
            for(n=0;n<p->gcin_count;++n)
            {
            i=p->gcin[n][0]+q;
            j=p->gcin[n][1];
            k=p->gcin[n][2];
            
            d->eddyv[IJK]=MIN(d->eddyv[IJK],1.0e-4);
            }
         pgc->start4V(p,d->eddyv,d->bc,24);
		}
        
}