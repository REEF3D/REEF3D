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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<iostream>

pic_f::pic_f(lexer *p, fdm *a, ghostcell *pgc, reini *ppreini, ioflow *ppflow, heat *pheat, concentration *pconc)
{
    static_assert(sizeof(pic_data)==NF*sizeof(double),"pic_data must be plain doubles");

    // the subdomain and the domain
    xloc_s = p->XN[0+marge];
    xloc_e = p->XN[p->knox+marge];
    yloc_s = p->YN[0+marge];
    yloc_e = p->YN[p->knoy+marge];
    zloc_s = p->ZN[0+marge];
    zloc_e = p->ZN[p->knoz+marge];

    xglo_s = pgc->globalmin(xloc_s);
    xglo_e = pgc->globalmax(xloc_e);
    yglo_s = pgc->globalmin(yloc_s);
    yglo_e = pgc->globalmax(yloc_e);
    zglo_s = pgc->globalmin(zloc_s);
    zglo_e = pgc->globalmax(zloc_e);

    fmarg = p->margin;
    npd = p->F83>=1 ? p->F83 : 1;
    npc = npd*npd*(p->j_dir==1 ? npd : 1);

    printtime = 0.0;
    printcount = 0;
    seeded = lost = lost_total = sent = sent_total = 0;
    added = removed = added_total = removed_total = 0;
    phitime = 0.0;

    // F 81 1, 2: the free surface from the particles, in place of the level set transport;
    // the density and viscosity follow it (the fluid update of the level set)
    surface = (p->F81==1 || p->F81==2);
    flipmode = (p->F81==2);
    relaxzones = surface && p->B90>0 && (p->B98==2 || p->B99==1 || p->B99==2);
    inflowsrc = surface && p->B60>0;
    preini = ppreini;
    pflow = ppflow;
    pupdate = nullptr;
    dfull = 1.0;
    uw = uo = km1 = nullptr;
    vw = vo = km2 = nullptr;
    ww = wo = km3 = nullptr;

    gcval_phi = 51;
    if(p->F50==2)
    gcval_phi = 52;
    if(p->F50==3)
    gcval_phi = 53;
    if(p->F50==4)
    gcval_phi = 54;

    if(surface)
    {
        if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
        pupdate = new fluid_update_fsf(p,a,pgc);

        if(p->F30>0 && p->H10==0 && p->W30==1 && p->F300==0 && p->W90==0)
        pupdate = new fluid_update_fsf_comp(p,a,pgc);

        if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==1)
        pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);

        if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==2)
        pupdate = new fluid_update_fsf_heat_Bouss(p,a,pgc,pheat);

        if(p->F30>0 && p->C10>0 && p->W90==0 && p->F300==0)
        pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);

        if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
        pupdate = new fluid_update_rheology(p);

        if(pupdate==nullptr)
        pupdate = new fluid_update_void();

        uw = new field1(p);
        vw = new field2(p);
        ww = new field3(p);
        km1 = new field1(p);
        km2 = new field2(p);
        km3 = new field3(p);
    }

    if(flipmode)
    {
        uo = new field1(p);
        vo = new field2(p);
        wo = new field3(p);
    }

    if(p->mpirank==0)
    {
    mkdir("./REEF3D_CFD_PIC",0777);

    logout.open("./REEF3D_CFD_PIC/REEF3D-CFD-PIC-Log.dat");
    logout<<"t \t particles \t left_domain \t sent_between_ranks \t added \t removed \t phi_time"<<endl;

    if(flipmode)
    cout<<"PIC (F 81 "<<p->F81<<"): the free surface and the momentum of the water from the particles (FLIP "<<p->F82<<"), "<<npc<<" per cell"<<endl;
    else if(surface)
    cout<<"PIC (F 81 "<<p->F81<<"): the free surface from the particles, "<<npc<<" per cell"<<endl;
    else
    cout<<"PIC (F 81 "<<p->F81<<"): passive particles, "<<npc<<" per cell"<<endl;
    }
}

pic_f::~pic_f()
{
    if(logout.is_open())
    logout.close();
}

void pic_f::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    seed(p,a,pgc);

    // B 60: the particles of the inflow in its ghost cells
    if(inflowsrc)
    inflow_seed(p,a,pgc);

    // the free surface of the start from the particles
    if(surface)
    {
    phi_ini(p,a,pgc);
    phi_build(p,a,pgc);
    water_velocity(p,a,pgc);
    }

    for(auto &q : P)
    velocity(p,a,q.x,q.y,q.z,q.u,q.v,q.w);

    // the first output (on a restart, I 45 1, too: the particles are seeded again)
    print(p,pgc);
}

void pic_f::start(lexer *p, fdm *a, ghostcell *pgc)
{
    // F 81 1, 2: the particles move with the velocity of the water
    if(surface)
    water_velocity(p,a,pgc);

    // F 81 2: the particles take the change of the water velocity
    if(flipmode)
    g2p(p,a);

    move(p,a,pgc);
    xchange(p,pgc);

    // F 78: the redistribution against the rows and gaps of the particles
    if(surface && p->F78>0 && p->count%p->F78==0)
    redistribute(p,a,pgc);

    if(surface)
    {
    phi_build(p,a,pgc);

    // the relaxation zones: the columns of the particles to the relaxed surface (pic_f_relax.cpp)
    if(relaxzones)
    zone_stretch(p,a);

    // the particles in the water (F 86) and, at every step, in the relaxation zones
    const bool interior = p->F86>0 && p->count%p->F86==0;

    if(interior || relaxzones)
    reseed(p,a,pgc,interior);

    if(inflowsrc)
    inflow_seed(p,a,pgc);
    }

    if(print_due(p))
    print(p,pgc);
}
