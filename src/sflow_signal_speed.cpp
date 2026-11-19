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

#include"sflow_signal_speed.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"

sflow_signal_speed::sflow_signal_speed(lexer* p) 
{

}

sflow_signal_speed::~sflow_signal_speed()
{
}

void sflow_signal_speed::signal_speed_update(lexer* p, ghostcell *pgc, fdm2D *b, 
                                        slice &Us, slice &Un, slice &Ve, slice &Vw, 
                                        slice &Ds, slice &Dn, slice &De, slice &Dw)
{
    const double g = fabs(p->W22);
    const int nj = p->jmax;
    const int *const wet = p->wet;
    
    // signal speed x-dir; the two wave speeds of a face are taken once (the same values as before)
    {
    const int *const flag = p->flagslice1;
    const double *const us = Us.V, *const un = Un.V, *const ds = Ds.V, *const dn = Dn.V;
    double *const Ss = b->Ss.V, *const Sn = b->Sn.V;
    const int ie = p->knox - p->ulast;
    
    for(int ii=0; ii<ie; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag[q]>0)
        {
        const double uL = us[q], uR = un[q];
        const double cL = sqrt(g*ds[q]), cR = sqrt(g*dn[q]);
        const int wL = wet[q], wR = wet[q+nj];
        
        if(wL==1 && wR==1)
        {
        const double US = 0.5*(uL+uR) + cL - cR;
        const double DS = 0.5*(cL + cR) + 0.25*(uL - uR);
        
        Ss[q] = MIN(uL - cL, US - DS);
        Sn[q] = MAX(uR + cR, US + DS);
        }
        
        else
        if(wL==0 && wR==1) // left dry
        {
        Ss[q] = uR - 2.0*cR;
        Sn[q] = uR + cR;
        }
        
        else
        if(wL==1 && wR==0) // right dry
        {
        Ss[q] = uL - cL;
        Sn[q] = uL + 2.0*cL;
        }
        
        else
        if(wL==0 && wR==0)
        {
        Ss[q] = 0.0;
        Sn[q] = 0.0;
        }
        }
    }
    }
    
    // signal speed y-dir
    if(p->j_dir==1)
    {
    const int *const flag = p->flagslice2;
    const double *const ve = Ve.V, *const vw = Vw.V, *const de = De.V, *const dw = Dw.V;
    double *const Se = b->Se.V, *const Sw = b->Sw.V;
    const int je = p->knoy - p->vlast;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+je; ++q)
        if(flag[q]>0)
        {
        const double uL = ve[q], uR = vw[q];
        const double cL = sqrt(g*de[q]), cR = sqrt(g*dw[q]);
        const int wL = wet[q], wR = wet[q+1];
        
        if(wL==1 && wR==1)
        {
        const double US = 0.5*(uL+uR) + cL - cR;
        const double DS = 0.5*(cL + cR) + 0.25*(uL - uR);
        
        Se[q] = MIN(uL - cL, US - DS);
        Sw[q] = MAX(uR + cR, US + DS);
        }
        
        else
        if(wL==0 && wR==1)
        {
        Se[q] = uR - 2.0*cR;
        Sw[q] = uR + cR;
        }
        
        else
        if(wL==1 && wR==0)
        {
        Se[q] = uL - cL;
        Sw[q] = uL + 2.0*cL;
        }
        
        else
        if(wL==0 && wR==0)
        {
        Se[q] = 0.0;
        Sw[q] = 0.0;
        }
        }
    }
    }
}
