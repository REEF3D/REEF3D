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

#include"sflow_forcing.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"nhflow_reinidisc_fsf.h"

sflow_forcing::sflow_forcing(lexer *p) : fx(p), fy(p), fz(p)
{
    forcing_flag=0;

    if(p->X10==2)
    forcing_flag=1;
}

sflow_forcing::~sflow_forcing()
{
}

void sflow_forcing::forcing(lexer *p, fdm2D *b, ghostcell *pgc, sixdof *p6dof, 
                             int iter, double alpha, slice &P, slice &Q, slice &w, slice &WL, bool finalize)
{
    if(forcing_flag==1)
    SLICELOOP4
    {
    fx(i,j) = 0.0;   
    fy(i,j) = 0.0;   
    }
    
    // 6DOF forcing
    p6dof->start_sflow(p,pgc,iter,b->fs,P,Q,w,b->fx,b->fy,b->fz,finalize);


    if(forcing_flag==1)
    {
    // add forcing term to RHS
    SLICELOOP1
    {
        P(i,j) += alpha*p->dt*fx(i,j);
        
        /*if(p->count<10)
        b->maxF = MAX(fabs(alpha*CPORNH*b->FX[IJK]), b->maxF);
        
        p->fbmax = MAX(fabs(alpha*CPORNH*b->FX[IJK]), p->fbmax);*/
    }
    
    SLICELOOP2
    {
        Q(i,j) += alpha*p->dt*fy(i,j);
        
        /*if(p->count<10)
        b->maxG = MAX(fabs(alpha*CPORNH*b->FY[IJK]), b->maxG);
        
        p->fbmax = MAX(fabs(alpha*CPORNH*b->FY[IJK]), p->fbmax);*/
    }
    
    SLICELOOP4
    {
        w(i,j) += alpha*p->dt*fz(i,j);
        
        /*if(p->count<10)
        b->maxG = MAX(fabs(alpha*CPORNH*b->FY[IJK]), b->maxG);
        
        p->fbmax = MAX(fabs(alpha*CPORNH*b->FY[IJK]), p->fbmax);*/
    }
    }
}

