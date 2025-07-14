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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"nhflow_reinidisc_fsf.h"

void nhflow_forcing::forcing_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(solid_flag==1)
    {
    if(p->mpirank==0)
    cout<<"Forcing ini "<<endl;
    
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
    
    pgc->start5V(p,p->ZSP,1);
    
    objects_create(p, pgc);
    ray_cast(p, d, pgc);
    reini_RK2(p, d, pgc, d->SOLID);
    
    SLICELOOP4
	d->depth(i,j) = p->wd - d->bed(i,j);
    }
    
    // --------------
    // DF
    LOOP
    p->DF[IJK]=1;
    
    if(solid_flag==1)
    LOOP
    if(d->SOLID[IJK]<0.0)
    p->DF[IJK]=-1;

    if(floating_flag==1)
    LOOP
    if(d->FB[IJK]<0.0)
    p->DF[IJK]=-1;
    
    pgc->startintV(p,p->DF,1);
    
    // -------------
    if(dlm_flag==1)
    objects_create(p, pgc);
    
    // DFSL slice
    pgc->gcsldf_update(p);
    pgc->solid_forcing_eta(p,d->WL);
    pgc->solid_forcing_eta(p,d->eta);
    pgc->solid_forcing_bed(p,d->bed);
}

void nhflow_forcing::reset(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(forcing_flag==1 || dlm_flag==1)
    {
    LOOP
    {
    FX[IJK] = 0.0;   
    FY[IJK] = 0.0;   
    FZ[IJK] = 0.0;   
    d->FHB[IJK] = 0.0;
    
    d->test[IJK] = 0.0;
    }

    SLICELOOP4
    fe(i,j) = 0.0;
    }
}
