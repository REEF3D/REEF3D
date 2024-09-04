/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "sediment_part.h"
#include "partres.h"
#include "lexer.h"
#include "ghostcell.h"
#include "fdm.h"
#include "vrans_f.h"
#include "reinitopo.h"
#include "ioflow.h"

/// @brief Updates the topography for the CFD solver
void sediment_part::update_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo* preto)
{
    prelax.start(p,pgc,&s);
    
    if(p->Q13==1)
    pst->update(p,*a,*pgc,PP);
    
    preto->start(p,a,pgc,a->topo);
    
    if(p->mpirank==0)
    cout<<"Topo: update grid..."<<endl;
    
    prelax.start(p,pgc,&s);
    pvrans->sed_update(p,a,pgc);
    pflow->gcio_update(p,a,pgc);
}
