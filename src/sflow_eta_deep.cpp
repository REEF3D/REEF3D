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

#include"sflow_eta.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sflow_eta::wetdrydeep(lexer* p, fdm2D* b, ghostcell* pgc, slice &eta, slice &P, slice &Q, slice &ws)
{
    //
    SLICELOOP4
    p->deep[IJ]=p->wet[IJ];
    
    SLICELOOP4
    {
    if(p->wet[Ip1J]==0 || p->wet[Ip2J]==0 || p->wet[Im1J]==0 || p->wet[Im2J]==0)
    p->deep[IJ]=0;
    
    
    if(p->j_dir==1)
    if(p->wet[IJp1]==0 || p->wet[IJp2]==0 || p->wet[IJp3]==0 || p->wet[IJm1]==0 || p->wet[IJm2]==0 || p->wet[IJm3]==0)
    p->deep[IJ]=0;
    
    if(p->j_dir==1)
    if(p->wet[Ip1Jp1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Im1Jm1]==0)
    p->deep[IJ]=0;
    }
    
    SLICELOOP4
    if(b->hp(i,j)<=10.0*wd_criterion)
    p->deep[IJ]=0;
    
    SLICELOOP1
    {
        if(p->deep[IJ]==1 && p->deep[Ip1J]==1)
        b->deep1(i,j)=1;
        
        if(p->deep[IJ]==0 || p->deep[Ip1J]==0)
        b->deep1(i,j)=0;
    }
    
    SLICELOOP2
    {
        if(p->deep[IJ]==1 && p->deep[IJp1]==1)
        b->deep2(i,j)=1;
        
        if(p->deep[IJ]==0 || p->deep[IJp1]==0)
        b->deep2(i,j)=0;
    }
    
    SLICELOOP4
    b->test(i,j) = p->deep[IJ];

    pgc->gcsl_start4Vint(p,p->deep,50);
    
}