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


#include"seastate_nhflow.h"
#include"seastate_f.h"
#include"seastate_vtp.h"
#include"fdm_seastate.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>

seastate_nhflow::seastate_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc) : seastate_coupling(p), Ub(p), Vb(p)
{
}

void seastate_nhflow::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->G1>0)
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  A 750 1 with NHFLOW: mesh refinement (G 1) is not supported"<<endl<<endl;

        pgc->final(true);
    }

    ini_wave(p,pgc,"NHFLOW");

    // surfbeat: the x- edge is open, its ghost cells are set here
    if(pwave->surfbeat())
    {
    p->open_xm = 1;

    ghostcells(p,d,d->UH,d->VH,d->WH);
    wl_ghostcells(p,d,d->WL);
    }

    average(p,d);
    forces(p,pgc,Ub,Vb);

    if(p->P10==1)
    pwave->printer()->print2D(p,pwave->e,pgc);
}

void seastate_nhflow::average(lexer *p, fdm_nhf *d)
{
    // depth-averaged velocity (sigma layers: sum DZN = 1), ghost cells included
    IMALOOP
    JMALOOP
    {
    double u=0.0, v=0.0;

        for(k=0; k<p->knoz; ++k)
        {
        u += d->U[IJK]*p->DZN[KP];
        v += d->V[IJK]*p->DZN[KP];
        }

    Ub(i,j) = u;
    Vb(i,j) = v*p->y_dir;
    }
}

void seastate_nhflow::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double dt;

    if(!due(p,dt))
    return;

    average(p,d);

    environment(p,d->eta,d->WL,Ub,Vb);

    wave_step(p,pgc,dt);

    forces(p,pgc,Ub,Vb);

    print(p,pgc);
}

void seastate_nhflow::u_source(lexer *p, fdm_nhf *d)
{
    const double r = ramp(p);

    LOOP
    WETDRY
    d->F[IJK] += r*Fw(i,j);
}

void seastate_nhflow::v_source(lexer *p, fdm_nhf *d)
{
    if(p->j_dir==0)
    return;

    const double r = ramp(p);

    LOOP
    WETDRY
    d->G[IJK] += r*Gw(i,j);
}

void seastate_nhflow::mass_source(lexer *p, fdm_nhf *d, slice &K)
{
    if(p->A751!=2)
    return;

    const double r = ramp(p);

    SLICELOOP4
    WETDRY
    K(i,j) -= r*divM(i,j);
}

void seastate_nhflow::ghostcells(lexer *p, fdm_nhf *d, double *UH, double *VH, double *WH)
{
    if(!pwave->surfbeat() || p->origin_i!=0)
    return;

    double etag, ug, vg;

    i = 0;

    for(j=0; j<p->knoy; ++j)
    {
        if(p->flagslice4[IJ]<0 || p->wet[IJ]==0)
        continue;

        double vbar = 0.0;
        for(k=0; k<p->knoz; ++k)
        vbar += d->V[IJK]*p->DZN[KP];

        if(!longwave_bc(p,j,d->WL(i,j),d->eta(i,j),vbar,etag,ug,vg))
        continue;

    const double hg = d->WL(i,j);

        for(k=0; k<p->knoz; ++k)
        {
        const double vk = d->V[IJK], wk = d->W[IJK];

        d->U[Im1JK]=d->U[Im2JK]=d->U[Im3JK]=ug;
        d->V[Im1JK]=d->V[Im2JK]=d->V[Im3JK]=vk;
        d->W[Im1JK]=d->W[Im2JK]=d->W[Im3JK]=wk;
        UH[Im1JK]=UH[Im2JK]=UH[Im3JK]=hg*ug;
        VH[Im1JK]=VH[Im2JK]=VH[Im3JK]=hg*vk;
        WH[Im1JK]=WH[Im2JK]=WH[Im3JK]=hg*wk;
        }
    }
}

void seastate_nhflow::wl_ghostcells(lexer *p, fdm_nhf *d, slice &WL)
{
    if(!pwave->surfbeat() || p->origin_i!=0)
    return;

    i = 0;

    for(j=0; j<p->knoy; ++j)
    {
        if(p->flagslice4[IJ]<0)
        continue;

        for(int q=1; q<=3; ++q)
        {
        WL(i-q,j) = WL(i,j);
        d->eta(i-q,j) = WL(i,j) - d->depth(i,j);
        }
    }
}
