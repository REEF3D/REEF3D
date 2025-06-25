/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"net_interface.h"
#include<sys/stat.h>

#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
    
sixdof_void::sixdof_void(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    mkdir("./REEF3D_CFD_6DOF",0777);
    
    pnetinter = new net_interface(p,pgc);
}

sixdof_void::~sixdof_void()
{
}

void sixdof_void::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_void::start_cfd(lexer* p, fdm* a, ghostcell* pgc, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{

    if(p->X310>0)
    {
        for (int i=0; i<p->mooring_count; i++)
        {
            pmooring[i]->start(p, pgc);
        }
    }
    
    if(p->X320>0)
    {
        pnetinter->netForces_cfd(p,a,pgc,alpha,Xne,Yne,Zne,Kne,Mne,Nne);
        
        NETLOOP
        {
        // Add to external forces
            Xext += Xne[n];
            Yext += Yne[n];
            Zext += Zne[n];
            Kext += Kne[n];
            Mext += Mne[n];
            Next += Nne[n];
        }
    }
    
    ++p->printcount_sixdof;
}


void sixdof_void::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}

void sixdof_void::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    
}

void sixdof_void::isource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_void::jsource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_void::ksource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_void::isource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_void::jsource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_void::ksource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_void::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void sixdof_void::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}
