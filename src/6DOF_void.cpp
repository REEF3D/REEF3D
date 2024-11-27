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
#include"vrans.h"
#include<sys/stat.h>

#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"
    
sixdof_void::sixdof_void(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    mkdir("./REEF3D_CFD_6DOF",0777);
}

sixdof_void::~sixdof_void()
{
}

void sixdof_void::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_void::start_cfd(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{

    if (p->X310 > 0)
    {
        for (int i=0; i<p->mooring_count; i++)
        {
            pmooring[i]->start(p, pgc);
        }
    }
    
    if (p->X320 > 0)
    {
        for (int ii = 0; ii < p->net_count; ii++)
        {
            pnet[ii]->start(p, a, pgc, 1.0, quatRotMat);
            pvrans->start(p, a, pgc, pnet[ii], ii);

            // Forces on rigid body
            pnet[ii]->netForces(p,Xne[ii],Yne[ii],Zne[ii],Kne[ii],Mne[ii],Nne[ii]);

            if( p->mpirank == 0)
            {
                cout<<"Xne"<< ii <<" : "<<Xne[ii]<<" Yne"<< ii <<" : "<<Yne[ii]<<" Zne"<< ii <<" : "<<Zne[ii]
                <<" Kne"<< ii <<" : "<<Kne[ii]<<" Mne"<< ii <<" : "<<Mne[ii]<<" Nne"<< ii <<" : "<<Nne[ii]<<endl;		
            }
        }
    }
    
    ++p->printcount_sixdof;
}


void sixdof_void::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}

void sixdof_void::start_sflow(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
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
