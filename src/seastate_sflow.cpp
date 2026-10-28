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


#include"patchBC_codes.h"
#include"seastate_sflow.h"
#include"seastate_remote.h"
#include"seastate_f.h"
#include"fdm_seastate.h"
#include"seastate_vtp.h"
#include"fdm2D.h"
#include"met_forcing.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>

seastate_sflow::seastate_sflow(lexer *p, fdm2D *b, ghostcell *pgc) : seastate_coupling(p)
{
}

void seastate_sflow::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the coupling acts on level 0: refined patches would run without the wave force and with law 5 for W 207 7
    if(p->G1>0)
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  A 750 with SFLOW: mesh refinement (G 1) is not supported (the patches get no"
            <<" wave force); for fine coastal cells use a stretched grid (DIVEMesh B 101 10, B 102 10, tools/stretched_grid.py)"<<endl<<endl;

        pgc->final(true);
    }

    // wave-dependent drag of the met forcing (W 207 7): the sea state of the coupling
    if(p->W207==7 && p->W200>0)
    {
        if(b->met!=nullptr)
        {
        b->met->waves(p,&Hsw,&Tpw,&b->WL);
        pmetw = b->met;
        }

        if(p->mpirank==0)
        cout<<(b->met!=nullptr ? "W 207 7: wave-dependent drag with the sea state of SEASTATE" : "!!! W 207 7: no met forcing on this grid, the drag stays law 5 !!!")<<endl;
    }

    // remote waves (A 750 2): SEASTATE is the other app of the MPI job
    if(p->A750==2)
    {
    seastate_remote::link = this;
    remote_ini(p,pgc,b->eta,b->WL,b->U,b->V);
    return;
    }

    ini_wave(p,pgc,"SFLOW");

    forces(p,pgc,b->U,b->V);

    if(p->P10>0)
    pwave->printer()->print2D(p,pwave->e,pgc);
}

void seastate_sflow::start(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double dt;

    if(remote)
    {
    remote_step(p,pgc,b->eta,b->WL,b->U,b->V,false);
    return;
    }

    if(!due(p,dt))
    return;

    environment(p,b->eta,b->WL,b->U,b->V);

    wave_step(p,pgc,dt);

    forces(p,pgc,b->U,b->V);

    print(p,pgc);
}

void seastate_sflow::u_source(lexer *p, fdm2D *b)
{
    const double r = ramp(p);

    SLICELOOP4
    WETDRY
    b->F(i,j) += r*Fw(i,j);
}

void seastate_sflow::v_source(lexer *p, fdm2D *b)
{
    const double r = ramp(p);

    SLICELOOP4
    WETDRY
    b->G(i,j) += r*Gw(i,j);
}

void seastate_sflow::mass_source(lexer *p, fdm2D *b, slice &K)
{
    if(p->A751!=2)
    return;

    const double r = ramp(p);

    SLICELOOP4
    WETDRY
    K(i,j) -= r*divM(i,j);
}

void seastate_sflow::ghostcells(lexer *p, fdm2D *b)
{
    if(!surfbeat() || p->origin_i!=0)
    return;

    double etag, ug, vg;

    GCSL4LOOP
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];

        if(p->gcbsl4[n][3]!=1 || patch_bc(p->gcbsl4[n][4]) || i!=0)
        continue;

        if(!longwave_bc(p,j,b->WL(i,j),b->eta(i,j),b->V(i,j),etag,ug,vg))
        continue;

        for(int q=1; q<=3; ++q)
        {
        b->eta(i-q,j) = etag;
        b->U(i-q,j) = ug;
        b->V(i-q,j) = vg;
        }
    }
}

void seastate_sflow::flux_bc(lexer *p, fdm2D *b, int ipol, slice &Fx)
{
    if(!surfbeat() || p->origin_i!=0)
    return;

    const double g = fabs(p->W22);

    GCSL4LOOP
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];

        if(p->gcbsl4[n][3]!=1 || patch_bc(p->gcbsl4[n][4]) || i!=0 || p->wet[IJ]==0)
        continue;

    const double eg = b->eta(i-1,j);
    const double wl = MAX(eg + b->depth(i-1,j), 0.0);
    const double ug = b->U(i-1,j);

        if(ipol==1)
        Fx(i-1,j) = wl*ug*ug + 0.5*g*eg*eg + g*eg*b->dfx(i-1,j);

        if(ipol==2)
        Fx(i-1,j) = wl*b->V(i-1,j)*ug;

        if(ipol==3)
        Fx(i-1,j) = wl*b->W(i-1,j)*ug;

        if(ipol==4)
        Fx(i-1,j) = wl*ug;
    }
}
