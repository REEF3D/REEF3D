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

#include"cfd_amr.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"rans_io.h"
#include"heat_print.h"
#include"concentration_io.h"
#include"cavitation.h"
#include<iostream>
#include<string>
#include<mpi.h>

//  The registry of the patch state: every field of a grid that the refinement moves between the
//  grids, with its kind (cells or the faces of a direction), its prolongation, its boundary values
//  and the parts of the step it takes part in (AF_*, cfd_amr.h).  The fills of the cells around
//  the patches, the restriction onto the coarser grids, the time levels of the subcycling, the
//  state of new patches (setup, regrid) and the boundary values of the grids are built from it.
//
//  The stage buffers of momentum_rk (the stage velocities, the level set of the stage input and
//  output) and the vectors of the composite projection are not state: the stage algorithm fills
//  and restricts them itself.  The density and viscosity follow from the level set
//  (margin_rovisc, rovisc_all); they are restricted, not filled.
//
//  A new field (heat T, concentration C, ...): one entry here with its roles, the field of the
//  patch's own object in get (patch_objects builds the object), its boundary code; the fills,
//  restriction, time levels and regrid transfer then take it.  AF_CONS marks a conserved quantity
//  (the heat or concentration of cfd_amr_scalar.cpp): its fluxes across the patch edges matched.
//
//  The order of the entries is the order of the boundary-value calls and of the time levels
//  (u, v, w first: their time levels are indexed by the direction).

void cfd_amr::register_fields()
{
    reg.clear();

    auto add = [&](const char *name, int type, int prolong, unsigned roles, int gc, int gcval,
                   std::function<field&(fdm*,turbulence*)> get)
    {
        cfd_amr_field e;
        e.name = name;
        e.type = type;
        e.prolong = prolong;
        e.roles = roles;
        e.gc = gc;
        e.gcval = gcval;
        e.get = get;
        reg.push_back(e);
    };

    const unsigned T = turb ? (AF_RSTEP|AF_RSTAGE|AF_TURBIN) : 0u;      // with a turbulence model
    const unsigned TB = turb ? AF_BCSTATE : 0u;

    // velocities (faces)
    const unsigned uvw = AF_STATE|AF_TIME|AF_RSTEP|AF_BCSTATE|AF_BCINI|AF_BCSTEP;
    add("u",0,1,uvw,1,10,[](fdm *a, turbulence*) -> field& { return a->u; });
    add("v",1,1,uvw,2,11,[](fdm *a, turbulence*) -> field& { return a->v; });
    add("w",2,1,uvw,3,12,[](fdm *a, turbulence*) -> field& { return a->w; });

    // level set, pressure
    add("phi",CELL,0,AF_STATE|AF_TIME|AF_STAGEIN|AF_RSTEP|AF_RSTAGE|AF_BCSTATE|AF_BCINI|AF_BCSTEP|AF_BCSTAGE,4,gcval_phi,
        [](fdm *a, turbulence*) -> field& { return a->phi; });
    add("press",CELL,0,AF_STATE|AF_TIME|AF_TMID|AF_STAGEIN|AF_RSTAGE|AF_BCSTATE|AF_BCINI|AF_BCSTAGE,4,40,
        [](fdm *a, turbulence*) -> field& { return a->press; });

    // eddy viscosity (zero without a model; its time levels and the stage fill always)
    add("eddyv",CELL,1,AF_STATE|AF_TIME|AF_STAGEIN|T|TB,4,24,[](fdm *a, turbulence*) -> field& { return a->eddyv; });

    // the RANS model: kin, eps (or omega); their boundary values by the model (turbulence::gcupdate)
    if(turb)
    {
        add("kin",CELL,1,AF_STATE|AF_TIME|T,0,0,[](fdm*, turbulence *t) -> field& { return dynamic_cast<rans_io*>(t)->kin; });
        add("eps",CELL,1,AF_STATE|AF_TIME|T,0,0,[](fdm*, turbulence *t) -> field& { return dynamic_cast<rans_io*>(t)->eps; });
        add("walld",CELL,0,AF_STATE,0,0,[](fdm *a, turbulence*) -> field& { return a->walld; });
    }

    // solids and topography (direct forcing): the level sets of the geometry
    if(solids)
    {
        add("solid",CELL,0,AF_STATE|AF_BCSTATE,5,150,[](fdm *a, turbulence*) -> field& { return a->solid; });
        add("topo",CELL,0,AF_STATE|AF_BCSTATE,5,150,[](fdm *a, turbulence*) -> field& { return a->topo; });
    }

    // density and viscosity: from the level set on every grid, restricted with it
    add("ro",CELL,0,AF_RSTEP|AF_RSTAGE|AF_BCINI|AF_BCSTEP|AF_BCSTAGE,4,1,[](fdm *a, turbulence*) -> field& { return a->ro; });
    add("visc",CELL,0,AF_RSTEP|AF_RSTAGE|AF_BCINI|AF_BCSTEP|AF_BCSTAGE,4,1,[](fdm *a, turbulence*) -> field& { return a->visc; });

    // heat T or concentration C (cfd_amr_scalar.cpp): conserved (flux registers), transported by
    // the refinement; the field of the heat / concentration object of the grid
    if(scal)
    {
        cfd_amr_field e;
        e.name = (scal==1) ? "T" : "C";
        e.type = CELL;
        e.prolong = 1;
        e.roles = AF_STATE|AF_TIME|AF_RSTEP|AF_BCSTATE|AF_BCINI|AF_BCSTEP|AF_CONS;
        e.gc = 4;
        e.gcval = 80;
        e.get = [](fdm *a, turbulence*) -> field& { return a->phi; };   // not used (pget, get0)
        e.pget = [this](cfd_amr_patch *c) -> field& { return scal==1 ? static_cast<heat_print*>(c->pheat)->T : static_cast<concentration_io*>(c->pconc)->C; };
        e.get0 = [this]() -> field& { return sfield(-1); };
        reg.push_back(e);
    }

    // cavitation (W 60): the vapour fraction of the cavitation object of the grid; filled around the
    // patches at every stage (the face densities of momentum and projection at the patch edges),
    // restricted, new patches from their parent; not flux-matched (cfd_amr_cav.cpp)
    if(cav)
    {
        cfd_amr_field e;
        e.name = "alpha_v";
        e.type = CELL;
        e.prolong = 1;
        e.roles = AF_STATE|AF_TIME|AF_STAGEIN|AF_RSTEP|AF_RSTAGE|AF_BCSTATE|AF_BCINI|AF_BCSTEP;
        e.gc = 4;
        e.gcval = 1;
        e.get = [](fdm *a, turbulence*) -> field& { return a->phi; };   // not used (pget, get0)
        e.pget = [](cfd_amr_patch *c) -> field& { return c->pcav->alpha; };
        e.get0 = [this]() -> field& { return cfield(-1); };
        reg.push_back(e);
    }

    // the time levels: u, v, w first (indexed by the direction)
    tl_n = 0;
    tl_phi = tl_press = tl_s = -1;
    for(cfd_amr_field &e : reg)
    if(e.roles & AF_TIME)
    {
        e.tl = tl_n++;
        if(e.type>=0 && e.tl!=e.type)
        {
            std::cout<<"CFD AMR: the registry needs u, v, w as the first time levels"<<std::endl;
            MPI_Abort(MPI_COMM_WORLD,-3111);
        }
        if(std::string(e.name)=="phi")
        tl_phi = e.tl;
        if(e.roles & AF_TMID)
        tl_press = e.tl;
        if(e.roles & AF_CONS)
        tl_s = e.tl;
    }

    // a field filled at a stage of G 7 1 is taken from the time levels of the coarser grid
    for(const cfd_amr_field &e : reg)
    if((e.roles & AF_STAGEIN) && !(e.roles & AF_TIME))
    {
        std::cout<<"CFD AMR: the registry field "<<e.name<<" is filled at a stage but has no time levels"<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,-3112);
    }
}

// the fills of the fields with the role
vector<cfd_amr::fspec> cfd_amr::specs(unsigned role)
{
    vector<fspec> fs;
    for(const cfd_amr_field &e : reg)
    if(e.roles & role)
    {
        const cfd_amr_field *pe = &e;
        fs.push_back({e.type,e.prolong,[this,pe](int g) -> field& { return rf(*pe,g); }});
    }
    return fs;
}

// the cell fields restricted with the role
vector<std::function<field&(int)>> cfd_amr::rcells(unsigned role)
{
    vector<std::function<field&(int)>> fs;
    for(const cfd_amr_field &e : reg)
    if((e.roles & role) && e.type==CELL)
    {
        const cfd_amr_field *pe = &e;
        fs.push_back([this,pe](int g) -> field& { return rf(*pe,g); });
    }
    return fs;
}

// the face fields restricted with the role, by direction (u, v, w)
vector<std::function<field&(int)>> cfd_amr::rfaces(unsigned role)
{
    vector<std::function<field&(int)>> fs(3);
    int nf = 0;
    for(const cfd_amr_field &e : reg)
    if((e.roles & role) && e.type>=0)
    {
        const cfd_amr_field *pe = &e;
        fs[e.type] = [this,pe](int g) -> field& { return rf(*pe,g); };
        ++nf;
    }
    if(nf==0)
    fs.clear();
    return fs;
}

// the boundary values of the fields with the role on grid g (the caller sets the scope of a patch)
void cfd_amr::bc_grid(int g, unsigned role)
{
    lexer *q = glex(g);
    for(const cfd_amr_field &e : reg)
    if((e.roles & role) && e.gc>0)
    {
        field &f = rf(e,g);
        if(e.gc==1) pgc0->start1(q,f,e.gcval);
        if(e.gc==2) pgc0->start2(q,f,e.gcval);
        if(e.gc==3) pgc0->start3(q,f,e.gcval);
        if(e.gc==4) pgc0->start4(q,f,e.gcval);
        if(e.gc==5) pgc0->start4a(q,f,e.gcval);
    }
}

// the kinds of the time levels
vector<int> cfd_amr::tl_types()
{
    vector<int> types;
    for(const cfd_amr_field &e : reg)
    if(e.roles & AF_TIME)
    types.push_back(e.type);
    return types;
}
