/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bedconc_VR.h"
#include"topo_relax.h"
#include"sediment_exnerdisc.h"
#include"sediment_fdm.h"

void sediment_exner::filter(lexer *p,ghostcell *pgc, slice &f, int outer_iter, int inner_iter)
{
    // Conservative, mask-aware smoothing of the bed velocity f (vz).
    //
    // Only sediment cells (DFBED>0, flagslice4>0) inside the S77 window take part.
    // Faces to structure/inactive cells, to cells outside the window and to the
    // domain boundary carry no flux (zero-gradient), so no bed change leaks into
    // cells where it is never applied.
    //
    // The face exchange  F = beta*min(A_i,A_n)*(h_n - h_i)  is symmetric, so
    // sum(A*f) is preserved on non-uniform grids. On a uniform grid this reduces
    // to the original  f = S102*h + (1-S102)*0.25*sum(h_n).
    //
    // predictor:  f = S(h)
    // corrector:  f += S(h - f), inner_iter times  (Van Cittert, sharpens the filter)

    slice4 h(p),dh(p),sdh(p),m(p);
    slice4 wxm(p),wxp(p),wym(p),wyp(p);

    const double beta = 0.25*(1.0-p->S102);
    const int ydir = (p->j_dir==1 && p->gknoy>1)?1:0;
    double Ai;

    // mask
    SLICEBASELOOP
    m(i,j)=0.0;

    SEDSLICELOOP
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    m(i,j)=1.0;

    pgc->gcsl_start4(p,m,1);

    // face weights
    SLICEBASELOOP
    {
    wxm(i,j)=wxp(i,j)=wym(i,j)=wyp(i,j)=0.0;

        if(m(i,j)>0.5)
        {
        Ai = p->DXN[IP]*p->DYN[JP];

        if(i+p->origin_i>0 && m(i-1,j)>0.5 && p->flagslice4[Im1J]>0 && p->DFBED[Im1J]>0)
        wxm(i,j) = beta*MIN(Ai,p->DXN[IM1]*p->DYN[JP])/Ai;

        if(i+p->origin_i<p->gknox-1 && m(i+1,j)>0.5 && p->flagslice4[Ip1J]>0 && p->DFBED[Ip1J]>0)
        wxp(i,j) = beta*MIN(Ai,p->DXN[IP1]*p->DYN[JP])/Ai;

            if(ydir==1)
            {
            if(j+p->origin_j>0 && m(i,j-1)>0.5 && p->flagslice4[IJm1]>0 && p->DFBED[IJm1]>0)
            wym(i,j) = beta*MIN(Ai,p->DXN[IP]*p->DYN[JM1])/Ai;

            if(j+p->origin_j<p->gknoy-1 && m(i,j+1)>0.5 && p->flagslice4[IJp1]>0 && p->DFBED[IJp1]>0)
            wyp(i,j) = beta*MIN(Ai,p->DXN[IP]*p->DYN[JP1])/Ai;
            }
        }
    }

	for(int qn=0;qn<outer_iter;++qn)
	{
		SLICEBASELOOP
		h(i,j) = f(i,j);

		pgc->gcsl_start4(p,h,1);

        // predictor
		SLICEBASELOOP
        if(m(i,j)>0.5)
		f(i,j) = h(i,j) + wxm(i,j)*(h(i-1,j)-h(i,j)) + wxp(i,j)*(h(i+1,j)-h(i,j))
                        + wym(i,j)*(h(i,j-1)-h(i,j)) + wyp(i,j)*(h(i,j+1)-h(i,j));

        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            SLICEBASELOOP
            dh(i,j) = h(i,j) - f(i,j);

            pgc->gcsl_start4(p,dh,1);

            SLICEBASELOOP
            {
            sdh(i,j) = 0.0;

            if(m(i,j)>0.5)
            sdh(i,j) = dh(i,j) + wxm(i,j)*(dh(i-1,j)-dh(i,j)) + wxp(i,j)*(dh(i+1,j)-dh(i,j))
                               + wym(i,j)*(dh(i,j-1)-dh(i,j)) + wyp(i,j)*(dh(i,j+1)-dh(i,j));
            }

            SLICEBASELOOP
            if(m(i,j)>0.5)
            f(i,j) += sdh(i,j);
		}
    }
}
