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
        count=0;
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];		

        uvel=uval[count]*ramp(p);
        vvel=vval[count]*ramp(p);
        wvel=wval[count]*ramp(p);
        
            U[Im1JK]=uvel;
            U[Im2JK]=uvel;
            U[Im3JK]=uvel;
            
            V[Im1JK]=vvel;
            V[Im2JK]=vvel;
            V[Im3JK]=vvel;
            
            W[Im1JK]=wvel;
            W[Im2JK]=wvel;
            W[Im3JK]=wvel;
            
        uvel=UHval[count]*ramp(p);
        vvel=VHval[count]*ramp(p);
        wvel=WHval[count]*ramp(p);
            
            
            UH[Im1JK]=uvel;
            UH[Im2JK]=uvel;
            UH[Im3JK]=uvel;
            
            VH[Im1JK]=vvel;
            VH[Im2JK]=vvel;
            VH[Im3JK]=vvel;
            
            WH[Im1JK]=wvel;
            WH[Im2JK]=wvel;
            WH[Im3JK]=wvel;
            
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
            
            d->EV[IJK]=MIN(d->EV[IJK],1.0e-4);
            }
        pgc->start4V(p,d->EV,24);
		}
}
