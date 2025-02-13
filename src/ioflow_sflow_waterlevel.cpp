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

#include"ioflow_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include"fdm2D.h"
#include"patchBC_interface.h"

void ioflow_f::waterlevel2D(lexer *p, fdm2D *b, ghostcell* pgc, slice &eta)
{
    // -------------------------------------
    // set fsf 
    double wsfout=p->phimean;
    double f=1.0;
    
    if(p->F62>1.0e-20)
    {
        if(p->F64==0)
        wsfout=p->F62;
        
        if(p->F64>0)
        {
        if(p->count<p->F64)
        f = 0.5*cos(PI + PI*double(p->count)/double(p->F64)) + 0.5;
        
        if(p->count>=p->F64)
        f = 1.0;
        
        wsfout = f*p->F62 + (1.0-f)*p->F60;
        }
    }
/*
        GCSL1LOOP
        {
        i = p->gcbsl1[n][0];
        j = p->gcbsl1[n][1];
                
            if(p->gcbsl1[n][4]==2)
            {
            //b->hx(i,j) = MAX(wsfout - b->bed(i,j),0.0);
            b->hx(i+1,j) = MAX(wsfout - b->bed(i,j),0.0);
            b->hx(i+2,j) = MAX(wsfout - b->bed(i,j),0.0);
            b->hx(i+3,j) = MAX(wsfout - b->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];
            
            if(p->gcbsl2[n][4]==2)
            {
            b->hy(i+1,j) = MAX(wsfout - b->bed(i,j),0.0);
            b->hy(i+2,j) = MAX(wsfout - b->bed(i,j),0.0);
            b->hy(i+3,j) = MAX(wsfout - b->bed(i,j),0.0);
            }
        }*/
        
        for(n=0;n<p->gcslout_count;n++)
        {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];
        
            if(p->wet[IJ]==1)
            {
            eta(i,j)   = wsfout-p->wd;
            eta(i+1,j) = wsfout-p->wd;
            eta(i+2,j) = wsfout-p->wd;
            eta(i+3,j) = wsfout-p->wd;

            b->hp(i,j)   = MAX(eta(i+1,j) + p->wd - b->bed(i,j),0.0);
            b->hp(i+1,j) = MAX(eta(i+1,j) + p->wd - b->bed(i,j),0.0);
            b->hp(i+2,j) = MAX(eta(i+2,j) + p->wd - b->bed(i,j),0.0);
            b->hp(i+3,j) = MAX(eta(i+3,j) + p->wd - b->bed(i,j),0.0);
            }
        }

    
    
    pBC->patchBC_waterlevel2D(p,b,pgc,eta);
}