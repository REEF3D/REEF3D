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

#include"fnpf_ini.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"reini.h"
#include"convection.h"
#include"ioflow.h"

fnpf_ini::fnpf_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_fsf_update(p,c,pgc), fnpf_bed_update(p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

fnpf_ini::~fnpf_ini()
{
}

void fnpf_ini::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, reini *preini, onephase *poneph)
{	
    
    //pflow->fi_relax(p,pgc,a->Fi,a->phi);
    //pflow->fifsf_relax(p,pgc,a->Fifsf);
    //pgc->start4(p,a->Fi,250);
    
    //pgc->gcsl_start4(p,a->Fifsf,50);
    pgc->gcsl_start4(p,c->eta,50);
}

void fnpf_ini::lsm_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow)
{
}

void fnpf_ini::velcalc(lexer *p, fdm_fnpf *c, ghostcell *pgc, field &f)
{
   
}
