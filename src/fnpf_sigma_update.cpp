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

#include"fnpf_sigma.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_fsf.h"

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e-5) // keep as is for wetting-drying

#define WLVLDRY (0.01*c->wd_criterion)


void fnpf_sigma::sigma_update(lexer *p, fdm_fnpf *c, ghostcell *pgc, fnpf_fsf *pf, slice &eta)
{
    // One fused pass per (i,j) column. All slice quantities are loaded once per
    // column; the per-cell arithmetic is kept in the original order so the
    // result is bit-identical to the previous five separate 3D sweeps.
    double *const __restrict sig   = p->sig;
    double *const __restrict sigx  = p->sigx;
    double *const __restrict sigy  = p->sigy;
    double *const __restrict sigxx = p->sigxx;
    double *const __restrict ZSN   = p->ZSN;
    double *const __restrict ZSP   = p->ZSP;
    const int knoz = p->knoz;
    
    ILOOP
    JLOOP
    {
        const double W   = WLVL;
        const double W2  = W*W;
        const double bx  = c->Bx(i,j),  by  = c->By(i,j);
        const double ex  = c->Ex(i,j),  ey  = c->Ey(i,j);
        const double bxx = c->Bxx(i,j), byy = c->Byy(i,j);
        const double exx = c->Exx(i,j), eyy = c->Eyy(i,j);
        
        const double bxW = bx/W, exW = ex/W;
        const double byW = by/W, eyW = ey/W;
        
        const double tbx = bxx - bx*bx/W;
        const double tex = exx - ex*ex/W;
        const double tby = byy - by*by/W;
        const double tey = eyy - ey*ey/W;
        const double sbex = bx + ex, bex = bx*ex;
        const double sbey = by + ey, bey = by*ey;
        
        k=0;
        const int n0 = FIJK;
        
        for(int kk=0; kk<=knoz; ++kk)
        {
            const int n = n0 + kk;
            const double s = sig[n];
            
            const double sgx = (1.0 - s)*bxW - s*exW;
            const double sgy = (1.0 - s)*byW - s*eyW;
            
            sigx[n] = sgx;
            sigy[n] = sgy;
            
            sigxx[n] = ((1.0 - s)/W)*tbx
                     - (s/W)*tex
                     - (sgx/W)*sbex
                     - ((1.0 - 2.0*s)/W2)*bex
                     + ((1.0 - s)/W)*tby
                     - (s/W)*tey
                     - (sgy/W)*sbey
                     - ((1.0 - 2.0*s)/W2)*bey;
        }
        
        // sigz
        if(p->flagslice4[IJ]>0)
        p->sigz[IJ] = 1.0/W;
        
        if(p->flagslice4[IJ]<0)
        p->sigz[IJ] = 1.0/(p->wd-p->bed[IJ]);
        
        // sig BC (was SLICELOOP4)
        if(p->flagslice4[IJ]>0)
        {
            const int nb = n0, nt = n0 + knoz;
            
            sigx[nb-1] = sigx[nb-2] = sigx[nb-3] = sigx[nb];
            sigx[nt+1] = sigx[nt+2] = sigx[nt+3] = sigx[nt];
            
            sigy[nb-1] = sigy[nb-2] = sigy[nb-3] = sigy[nb];
            sigy[nt+1] = sigy[nt+2] = sigy[nt+3] = sigy[nt];
            
            sigxx[nb-1] = sigxx[nb-2] = sigxx[nb-3] = sigxx[nb];
            sigxx[nt+1] = sigxx[nt+2] = sigxx[nt+3] = sigxx[nt];
        }
        
        // ZSN (was FLOOP) and ZSP (was LOOP)
        const double wl = c->WL(i,j), bd = c->bed(i,j);
        
        for(k=0; k<=knoz; ++k)
        {
            // FLOOP already implies flag7>0, so the old FSCHECK branch was dead
            if(p->flag7[FIJK]>0)
            ZSN[FIJK] = p->ZN[KP]*wl + bd;
        }
        
        for(k=0; k<knoz; ++k)
        if(p->flag4[IJK]>0)
        ZSP[IJK]  = p->ZP[KP]*wl + bd;
    }
}
