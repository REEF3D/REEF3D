/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans_v.h"
#include"vrans_f.h"
#include"vrans_veg.h"

void ioflow_f::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    gcio_update(p,a,pgc);
    
    
    if(p->B269==0)
	pvrans = new vrans_v(p,a,pgc);
	
	if(p->B269==1 || p->S10==2)
	pvrans = new vrans_f(p,a,pgc);
    
    if(p->B269==2)
	pvrans = new vrans_veg(p,a,pgc);
}

void ioflow_f::ini_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void ioflow_f::ini2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
}

void ioflow_f::full_initialize2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

