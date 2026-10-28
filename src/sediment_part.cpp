/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"lexer.h"
#include"bedshear.h"
#include"momentum_forcing.h"
#include"CPM.h"
#include"suspended.h"
#include"sediment_exner.h"
#include"ckpt_pack.h"
// the classes deleted in the destructor (deleting a pointer to an incomplete type ran no destructor)
#include"sediment_fdm.h"
#include"bedload.h"
#include"bedconc.h"
#include"sandslide.h"
#include"topo_relax.h"
#include"vrans.h"
#include"bedslope.h"
#include"bedshear_reduction.h"
#include"topo.h"
#include"diffusion.h"
#include"convection.h"
#include"bedload_direction.h"
#include<sys/stat.h>

sediment_part::sediment_part(lexer *p, fdm *a, ghostcell *pgc, turbulence *ppturb, patchBC_interface *ppBC) : por(p), d50(p)
{
    pBC = ppBC;
    pturb = ppturb;

    sediment_logic(p,pgc);
    
    // two-way coupling: momentum source of the parcels on the fluid
    if(p->Q50==1 && p->Q11==2)
    momentum_forcing::psed = this;

    // Create Folder
    if(p->mpirank==0 && p->Q180>0 && (p->Q181>0||p->Q182>0))
        mkdir("./REEF3D_CFD_CPM_Particle",0777);
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
    delete pbeddir;
    delete pslope;
    // pBC and pturb belong to the caller (constructor arguments): they were deleted here as well
}

// hybrid suspension (Q 58 3, see CPM_suspension.cpp): erosion concentration from the parcel bed,
// concentration step (suspended_IM1), net exchange with the bed for the parcels
void sediment_part::start_susp(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
    if(psusp==nullptr || p->Q44==1)
    return;
    
    if(!((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47)))
    return;
    
    pst->susp_cbe(p,a,pgc,s);
    psusp->start(a,p,psuspdisc,psuspdiff,psolv,pgc,pflow,s);
    pst->susp_flux(p,a,pgc,s);
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
    
    // two-way coupling once the parcels move (S 41, S 43/S 45: the fluid spins up first)
    if(!sedcalc)
    pst->coupling_update(p,a,pgc,s);
    else
    pst->coupling_reset(p,pgc);
    
    print_particles(p,s);
    
}

// checkpoints (P 48, I 45): the scalars of the bed update and of the wave-averaged bed shear (S 16 9);
// the fields and slices (the CPM fields, the bed) are in the registry of the grid, the parcels in the
// rank blob (CPM::ckpt_rank_state)
void sediment_part::ckpt_state(lexer *p, std::vector<double> &v)
{
    ckpt_put w{v};
    w.d(1.0);   // version
    std::vector<double> ex;
    static_cast<sediment_exner*>(ptopo)->ckpt_state(ex);
    w.vd(ex);
    w.d(*pbedshear->ckpt_tw_old()); w.d(*pbedshear->ckpt_tw_dt()); w.i(*pbedshear->ckpt_tw_count());
}

bool sediment_part::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
{
    ckpt_get r{v};
    double version=0.0;
    r.d(version);
    if(!r.ok || int(version)!=1)
    return false;
    std::vector<double> ex;
    r.vd(ex);
    size_t k=0;
    static_cast<sediment_exner*>(ptopo)->ckpt_restore(ex,k);
    r.d(*pbedshear->ckpt_tw_old()); r.d(*pbedshear->ckpt_tw_dt()); r.i(*pbedshear->ckpt_tw_count());
    return r.ok && k==ex.size();
}

void sediment_part::ckpt_rank_state(lexer *p, std::vector<double> &v)
{
    pst->ckpt_rank_state(p,v);
}

bool sediment_part::ckpt_rank_restore(lexer *p, ghostcell *pgc, const std::vector<std::vector<double>> &v, bool same)
{
    return pst->ckpt_rank_restore(p,pgc,v,same);
}

// G 46: the CPM load map (parcels per cell) for the weighted partition of the restart
void sediment_part::loadmap(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->A10==6)
    pst->loadmap_write(p,a,pgc);
}
