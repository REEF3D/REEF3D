/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"lexer.h"
#include"partres.h"
#include"sediment_fdm.h"
#include"sandslide_f.h"
#include"sandslide_f2.h"
#include"sandslide_f3.h"
#include"sandslide_pde.h"
#include"sandslide_v.h"
#include"vrans_f.h"
#include"vrans_v.h"
#include"bedslope.h"
#include"reduction_parker.h"
#include"reduction_deyemp.h"
#include"reduction_deyana.h"
#include"reduction_FD.h"
#include"reduction_void.h"
#include"sediment_exner.h"
#include"topo_relax.h"
#include"bedshear.h"

void sediment_part::sediment_logic(lexer *p, ghostcell *pgc)
{
    pst = new partres(p,pgc);

    s = new sediment_fdm(p);

    if(p->S90==1)
        pslide = new sandslide_f(p);
    else if(p->S90==2)
        pslide = new sandslide_f2(p);
    else if(p->S90==3)
        pslide = new sandslide_f3(p);
    else if(p->S90==4)
        pslide = new sandslide_pde(p);
    else
        pslide = new sandslide_v(p);

    if(p->S10==2 && p->A10==6)
        pvrans = new vrans_f(p,pgc);
    else
        pvrans = new vrans_v(p,pgc);

    pslope = new bedslope(p);

    if(p->S80==1)
        preduce = new reduction_parker(p);
    else if(p->S80==2)
        preduce = new reduction_deyemp(p);
    else if(p->S80==3)
        preduce = new reduction_deyana(p);
    else if(p->S80==4)
        preduce = new reduction_FD(p);
    else
        preduce = new reduction_void(p);

    ptopo = new sediment_exner(p,pgc);

    prelax = new topo_relax(p);

    pbedshear = new bedshear(p,pturb);
}
