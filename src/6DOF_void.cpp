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
    
    alpha[0] = 8.0/15.0;
    alpha[1] = 2.0/15.0;
    alpha[2] = 2.0/6.0;
    
    gamma[0] = 8.0/15.0;
    gamma[1] = 5.0/12.0;
    gamma[2] = 3.0/4.0;
    
    zeta[0] = 0.0;
    zeta[1] = -17.0/60.0;
    zeta[2] = -5.0/12.0;
    
    if(((p->N40==3 || p->N40==13 || p->N40==23 || p->N40==33) && p->A10==6) || (p->A510==3 && p->A10==5) || (p->A210==3 && p->A10==2)) 
    {
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
    
    gamma[0] = 0.0;
    gamma[1] = 0.0;
    gamma[2] = 0.0;
    
    zeta[0] = 0.0;
    zeta[1] = 0.0;
    zeta[2] = 0.0;
    }
    
    if(((p->N40==2 || p->N40==12 || p->N40==22) && p->A10==6) || (p->A510==2 && p->A10==5) || (p->A210==2 && p->A10==2)) 
    {
    alpha[0] = 1.0;
    alpha[1] = 0.5;
    
    gamma[0] = 0.0;
    gamma[1] = 0.0;
    gamma[2] = 0.0;
    
    zeta[0] = 0.0;
    zeta[1] = 0.0;
    zeta[2] = 0.0;
    }
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
        pnetinter->netForces_cfd(p,a,pgc,alpha[iter],quatRotMat,Xne,Yne,Zne,Kne,Mne,Nne);
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
