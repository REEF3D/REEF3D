/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fnpf_sg_fsf_update.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"onephase.h"
#include"slice.h"

fnpf_sg_fsf_update::fnpf_sg_fsf_update(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
    gcval_u = 10;
    gcval_v = 11;
    gcval_w = 12;
}

fnpf_sg_fsf_update::~fnpf_sg_fsf_update()
{
    
}

void fnpf_sg_fsf_update::fsfupdate(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, onephase *poneph, slice &eta)
{
}

void fnpf_sg_fsf_update::etaloc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void fnpf_sg_fsf_update::etaloc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // find k location for eta
    SLICELOOP4
    c->etaloc(i,j) = p->knoz;
    
    pgc->gcsl_start4int(p,c->etaloc,50);
}

void fnpf_sg_fsf_update::fsfbc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Fifsf, double *Fi)
{
    FFILOOP4
    {
        //Fi[FIJK]   = Fifsf(i,j);
        Fi[FIJKp1] = Fifsf(i,j);
        Fi[FIJKp2] = Fifsf(i,j);  
        Fi[FIJKp3] = Fifsf(i,j);
    }
}

void fnpf_sg_fsf_update::fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Fifsf, field &Fi)
{
}

void fnpf_sg_fsf_update::fsfepol(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, field &Fi)
{
}

void fnpf_sg_fsf_update::velcalc(lexer *p, fdm_fnpf *c, ghostcell *pgc, field &f)
{
}

void fnpf_sg_fsf_update::velcalc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{
    
    LOOP
    c->Fi4(i,j,k) = 0.5*(c->Fi[FIJK]+c->Fi[FIJKp1]);
    
    pgc->start4(p,c->Fi4,250);
    
    LOOP
    {
    c->u(i,j,k) = (c->Fi4(i+1,j,k)-c->Fi4(i,j,k))/(p->DXP[IP]) + p->sigx[FIJK]*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZP[KP]+p->DZP[KM1]));
    }
    
    LOOP
    {
	c->v(i,j,k) = ((c->Fi4(i,j+1,k)-c->Fi4(i,j,k))/(p->DYP[JP]))+ p->sigy[FIJK]*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZP[KP]+p->DZP[KM1]));
    }
    
    LOOP
    {
    c->w(i,j,k) = ((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZP[KP]+p->DZP[KM1]))*p->sigz[FIJK];
    }
    
    pgc->start4(p,c->u,1);
	pgc->start4(p,c->v,1);
	pgc->start4(p,c->w,1);
    
    c->u.ggcpol(p);
	c->v.ggcpol(p);
	c->w.ggcpol(p);
}
