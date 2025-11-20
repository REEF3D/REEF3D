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
#include"bedshear.h"
#include<sys/stat.h>

sediment_part::sediment_part(lexer *p, fdm *a, ghostcell *pgc, turbulence *ppturb, patchBC_interface *ppBC) : por(p), d50(p)
{
    pBC = ppBC;
    pturb = ppturb;

    sediment_logic(p,pgc);

    // Create Folder
    if(p->mpirank==0 && p->Q180>0 && (p->Q181>0||p->Q182>0))
        mkdir("./REEF3D_CFD_SedPart",0777);
}

sediment_part::~sediment_part()
{
    delete pst;
    delete s;
    delete pbed;
    delete pvrans;
    delete pcbed;
    delete pslide;
    delete prelax;
    delete preduce;
    delete ptopo;
    delete psusp;
    delete psuspdiff;
    delete psuspdisc;
    delete pbedshear;
    delete pBC;
    delete pbeddir;
    delete pslope;
    delete pturb;
}

void sediment_part::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo *preto, solver *psolv)
{
    bool sedcalc = true;

    if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
    {
        sediment_algorithm_cfd(p,a,pgc,pflow,preto);

        sedcalc = false;
    }

    if(sedcalc)
    {
        fill_bedk(p,a,pgc);
        waterlevel(p,a,pgc);
        pbedshear->taubed(p,a,pgc,s);
    }
}
