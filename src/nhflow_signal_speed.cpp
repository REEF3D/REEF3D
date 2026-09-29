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
FITNESS FOR A PARTICULAR PURPOSE. d->Se[IJK]e the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_signal_speed.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_signal_speed::nhflow_signal_speed(lexer* p) 
{

}

nhflow_signal_speed::~nhflow_signal_speed()
{
}

void nhflow_signal_speed::signal_speed_update(lexer* p, ghostcell *pgc, fdm_nhf *d, 
                                        double *Us, double *Un, double *Ve, double *Vw, 
                                        slice &Ds,slice &Dn, slice &De, slice &Dw)
{
    // Same expressions as before; the celerities sqrt(9.81*D) depend only on
    // the column and are evaluated once per (i,j) instead of 4-6 times per cell.
    
    // signal speed x-dir  (ULOOP)
    for(i=0; i<p->knox-p->ulast; ++i)
    for(j=0; j<p->knoy; ++j)
    {
    const double cs = sqrt(9.81*Ds(i,j));
    const double cn = sqrt(9.81*Dn(i,j));
    const int wP = p->wet[IJ];
    const int wN = p->wet[Ip1J];
    
        for(k=0; k<p->knoz; ++k)
        if(p->flag1[IJK]>0)
        {
        const double us = Us[IJK];
        const double un = Un[IJK];
        const double usx = 0.5*(us+un) + cs - cn;
        const double dsx = 0.5*(cs + cn) + 0.25*(us - un);
        const int dfP = p->DF[IJK];
        const int dfN = p->DF[Ip1JK];

        if((wP==1 && wN==1) && (dfP==1 && dfN==1))
        {
        d->Ss[IJK] = MIN(us - cs, usx - dsx);
        d->Sn[IJK] = MAX(un + cn, usx + dsx);
        d->SSx[IJK] = usx;
        }
        
        else
        if((wP==0 && wN==1) || (dfP<0 && dfN==1))  // left dry
        {
        d->Ss[IJK] = un - 2.0*cn;
        d->Sn[IJK] = un +     cn;
        d->SSx[IJK] = d->Ss[IJK];
        }
        
        else
        if((wP==1 && wN==0) || (dfP==1 && dfN<0)) // right dry
        {
        d->Ss[IJK] = us - cs;
        d->Sn[IJK] = us + cs;
        d->SSx[IJK] = d->Sn[IJK];
        }
        
        else
        if((wP==0 && wN==0)  || (dfP<0 && dfN<0))
        {
        d->Ss[IJK] = 0.0;
        d->Sn[IJK] = 0.0;
        d->SSx[IJK] = 0.0;
        }
        }
    }
    
    // signal speed y-dir  (VLOOP)
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=0; j<p->knoy-p->vlast; ++j)
    {
    const double ce = sqrt(9.81*De(i,j));
    const double cw = sqrt(9.81*Dw(i,j));
    const int wP = p->wet[IJ];
    const int wN = p->wet[IJp1];
    
        for(k=0; k<p->knoz; ++k)
        if(p->flag2[IJK]>0)
        {
        const double ve = Ve[IJK];
        const double vw = Vw[IJK];
        const double usy = 0.5*(ve+vw) + ce - cw;
        const double dsy = 0.5*(ce + cw) + 0.25*(ve - vw);
        const int dfP = p->DF[IJK];
        const int dfN = p->DF[IJp1K];
        
        if((wP==1 && wN==1) && (dfP==1 && dfN==1))
        {
        d->Se[IJK] = MIN(ve - ce, usy - dsy);
        d->Sw[IJK] = MAX(vw + cw, usy + dsy);
        d->SSy[IJK] = usy;
        }

        else
        if((wP==0 && wN==1) || (dfP<0 && dfN==1))
        {
        d->Se[IJK] = vw - 2.0*cw;
        d->Sw[IJK] = vw +     cw;
        d->SSy[IJK] = d->Se[IJK];
        }
        
        else
        if((wP==1 && wN==0) || (dfP==1 && dfN<0))
        {
        d->Se[IJK] = ve -     ce;
        d->Sw[IJK] = ve + 2.0*ce;
        d->SSy[IJK] = d->Sw[IJK];
        }
        
        else
        if((wP==0 && wN==0) || (dfP<0 && dfN<0))
        {
        d->Se[IJK] = 0.0;
        d->Sw[IJK] = 0.0;
        d->SSy[IJK] = 0.0;
        }
        }
    }
}
