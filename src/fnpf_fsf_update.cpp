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

#include"fnpf_fsf_update.h"#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"slice.h"

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e20)

fnpf_fsf_update::fnpf_fsf_update(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
    gcval_u = 10;
    gcval_v = 11;
    gcval_w = 12;
}

fnpf_fsf_update::~fnpf_fsf_update()
{
    
}

void fnpf_fsf_update::fsfupdate(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, slice &eta)
{
}

void fnpf_fsf_update::etaloc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void fnpf_fsf_update::etaloc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // find k location for eta
    SLICELOOP4
    c->etaloc(i,j) = p->knoz;
    
    pgc->gcsl_start4int(p,c->etaloc,50);
}

void fnpf_fsf_update::fsfbc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Fifsf, double *Fi)
{
    FFILOOP4
    {
        Fi[FIJK]   = Fifsf(i,j);
        Fi[FIJKp1] = Fifsf(i,j);
        Fi[FIJKp2] = Fifsf(i,j);  
        Fi[FIJKp3] = Fifsf(i,j);
    }
}

void fnpf_fsf_update::fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Fifsf, field &Fi)
{
}

void fnpf_fsf_update::fsfepol(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, field &Fi)
{
}

void fnpf_fsf_update::velcalc(lexer *p, fdm_fnpf *c, ghostcell *pgc, field &f)
{
}

