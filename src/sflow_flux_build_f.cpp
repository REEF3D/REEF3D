/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sflow_flux_build_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"slice.h"
#include"patchBC_interface.h"

sflow_flux_build_f::sflow_flux_build_f(lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pBC = ppBC;
}

sflow_flux_build_f::~sflow_flux_build_f()
{
}

void sflow_flux_build_f::start_U(lexer* p, fdm2D *b, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    b->Fs(i,j) = b->UHs(i,j)*b->Us(i,j)
            + 0.5*fabs(p->W22)*b->ETAs(i,j)*b->ETAs(i,j)
            + fabs(p->W22)*b->ETAs(i,j)*b->dfx(i,j);
    
    b->Fn(i,j) = b->UHn(i,j)*b->Un(i,j)
            + 0.5*fabs(p->W22)*b->ETAn(i,j)*b->ETAn(i,j)
            + fabs(p->W22)*b->ETAn(i,j)*b->dfx(i,j);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    b->Fe(i,j) = b->Ve(i,j)*b->UHe(i,j);
    
    b->Fw(i,j) = b->Vw(i,j)*b->UHw(i,j);
    }
}

void sflow_flux_build_f::start_V(lexer* p, fdm2D *b, ghostcell *pgc)
{
    if(p->j_dir==1)
    {
    // flux x-dir
    ULOOP
    {
    b->Fs(i,j) = b->Us(i,j)*b->VHs(i,j);
    
    b->Fn(i,j) = b->Un(i,j)*b->VHn(i,j);
    }
    
    // flux y-dir
    VLOOP
    {
    b->Fe(i,j) = b->VHe(i,j)*b->Ve(i,j) 
            + 0.5*fabs(p->W22)*b->ETAe(i,j)*b->ETAe(i,j)
            + fabs(p->W22)*b->ETAe(i,j)*b->dfy(i,j);
    
    b->Fw(i,j) = b->VHw(i,j)*b->Vw(i,j) 
            + 0.5*fabs(p->W22)*b->ETAw(i,j)*b->ETAw(i,j) 
            + fabs(p->W22)*b->ETAw(i,j)*b->dfy(i,j);
    }
    }
}

void sflow_flux_build_f::start_W(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    b->Fs(i,j) = b->Us(i,j)*b->WHs(i,j);
    
    b->Fn(i,j) = b->Un(i,j)*b->WHn(i,j);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    b->Fe(i,j) = b->Ve(i,j)*b->WHe(i,j);
    
    b->Fw(i,j) = b->Vw(i,j)*b->WHw(i,j);
    }
}

void sflow_flux_build_f::start_E(lexer* p, fdm2D *b, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    b->Fs(i,j) = b->UHs(i,j);
    
    b->Fn(i,j) = b->UHn(i,j);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    b->Fe(i,j) = b->VHe(i,j);
    
    b->Fw(i,j) = b->VHw(i,j);
    }
}





