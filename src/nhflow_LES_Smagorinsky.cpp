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

#include"nhflow_LES_Smagorinsky.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_strain.h"
#include"solver.h"
#include"nhflow_diffusion.h"
#include"ioflow.h"
#include"nhflow_scalar_convection.h"

nhflow_LES_Smagorinsky::nhflow_LES_Smagorinsky(lexer* p, fdm_nhf* d, ghostcell *pgc) : nhflow_les_io(p,d)
{
    c_sgs=0.2;
}

nhflow_LES_Smagorinsky::~nhflow_LES_Smagorinsky()
{
}

void nhflow_LES_Smagorinsky::start(lexer* p, fdm_nhf* d, ghostcell* pgc, nhflow_scalar_convection* pconvec, nhflow_diffusion* pdiff,solver* psolv, ioflow* pflow, vrans *pvrans)
{
	LOOP
    d->EV[IJK] = pow(c_sgs,2.0) * pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*d->WL(i,j),2.0/3.0) * strainterm(p,d);
    
    
    double s11,s22,s33,s12,s13,s23;
    
    LOOP
    {
        s11 = dudx(d->U);
        s22 = dvdy(d->V);
        s33 = dwdz(d->W);
        s12 = (dudy(d->U) + dvdx(d->V));
        s13 = (dudz(d->U) + dwdx(d->W));
        s23 = (dvdz(d->V) + dwdy(d->W));
        
    d->test[IJK]=s11;
    }

    pgc->start24V(p,d->EV,24);
}

void nhflow_LES_Smagorinsky::ktimesave(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
}

void nhflow_LES_Smagorinsky::etimesave(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
}

void nhflow_LES_Smagorinsky::kinupdate(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
}

void nhflow_LES_Smagorinsky::timesource(lexer* p, fdm_nhf* d, double *FN)
{
}

void nhflow_LES_Smagorinsky::clearrhs(lexer* p, fdm_nhf *d)
{
}
