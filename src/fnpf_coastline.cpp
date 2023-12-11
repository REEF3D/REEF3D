/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"fnpf_coastline.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice.h"
#include"sliceint.h"

fnpf_coastline::fnpf_coastline(lexer* p) :  ddweno_f_nug(p), frk1(p),frk2(p),L(p),dt(p),wet_n(p)
{
    time_preproc(p); 
}

fnpf_coastline::~fnpf_coastline()
{
}

void fnpf_coastline::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &coastline, int *wet, sliceint &wet_n)
{
    if(p->count==0)
    {
        SLICELOOP4
        {

            coastline(i,j)=1.0;
            
            if(p->wd - c->bed(i,j) < c->wd_criterion)
            coastline(i,j)=-1.0;
   
        }
        reini(p,pgc,coastline);
    }
}




