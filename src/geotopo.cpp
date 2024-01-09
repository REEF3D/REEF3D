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

#include"geotopo.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinitopo.h"
#include"ioflow.h"

geotopo::geotopo(lexer* p, fdm *a, ghostcell* pgc)
{
}

geotopo::~geotopo()
{
}

void geotopo::start(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow, reinitopo* preto, vrans* pvrans)
{
    dat(p,a,pgc);
    
    preto->start(p,a,pgc,a->topo);
    
    if(p->G3==0)
    {
        if(p->S10!=2)
        pgc->topo_update(p,a);
        
        pflow->gcio_update(p,a,pgc);
    }
    
    if(p->S10==2)
    pflow->vrans_sed_update(p,a,pgc,pvrans);
}



