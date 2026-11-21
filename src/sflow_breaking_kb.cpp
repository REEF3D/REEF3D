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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sflow_breaking_kb.h"
#include"breaking_crest3d.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<vector>
#include<sys/stat.h>
#include<cstdlib>

sflow_breaking_kb::sflow_breaking_kb(lexer *p) : brkage(p), brkage_n(p), brkk(p), brkk_n(p), brkr(p), brkr_n(p), brkkseed(p),
                                                 us(p), vs(p), vbs(p), vbk(p), vbp(p), brkseed(p), brkflag(p), brkzone(p), bd(p)
{
    SLICELOOP4
    {
    brkage(i,j)=brkage_n(i,j)=brkk(i,j)=brkk_n(i,j)=brkr(i,j)=brkr_n(i,j)=brkkseed(i,j)=0.0;
    us(i,j)=vs(i,j)=vbs(i,j)=vbk(i,j)=vbp(i,j)=0.0;
    brkseed(i,j)=brkflag(i,j)=brkzone(i,j)=bd(i,j)=0;
    }
}

sflow_breaking_kb::~sflow_breaking_kb()
{
    delete bt;
}

void sflow_breaking_kb::zones(lexer *p, ioflow *pflow)
{
    if(zones_built)
    return;

    // boxes xs, xe, ys, ye of the wave generation and the numerical beach
    std::vector<double> bx;
    if(pflow!=nullptr)
    pflow->relax_boxes(p,bx);

    SLICELOOP4
    {
    brkzone(i,j) = 0;

        for(size_t q=0; q+3<bx.size(); q+=4)
        if(p->XP[IP]>=bx[q] && p->XP[IP]<=bx[q+1] && (p->j_dir==0 || (p->YP[JP]>=bx[q+2] && p->YP[JP]<=bx[q+3])))
        brkzone(i,j) = 1;
    }

    // the open boundaries of the wave generation (B 98 3/4) and the outflow: deta/dt there is the forcing, the
    // viscosity next to the boundary values made the explicit SFLOW unstable; no breaking within 5 cells
    for(int q=0; q<2; ++q)
    {
    const int nb = (q==0) ? p->gcslin_count : p->gcslout_count;

        for(int m=0; m<nb; ++m)
        {
        const int bi = (q==0) ? p->gcslin[m][0] : p->gcslout[m][0];
        const int bj = (q==0) ? p->gcslin[m][1] : p->gcslout[m][1];

        SLICELOOP4
        if(abs(i-bi)<=5 && abs(j-bj)<=5)
        brkzone(i,j) = 1;
        }
    }

    zones_built = true;

    int nz = 0;
    SLICELOOP4
    nz += brkzone(i,j);

    if(p->mpirank==0)
    cout<<"SFLOW A 246 4: breaking with crest tracking and the eddy viscosity of Kennedy et al. (2000), "<<nz<<" cells without breaking (relaxation zones, open boundaries) on rank 0"<<endl;
}

void sflow_breaking_kb::start(lexer *p, fdm2D *b, ghostcell *pgc, slice &eta, slice &eta_n, double alpha)
{
    // crest tracking once per time step (at the first call of the step: eta of the stage closest to t^n),
    // the viscosity at every call (deta/dt of the stage)
    tracking(p,b,pgc,eta);
    viscosity(p,b,pgc,eta,eta_n,alpha);
}

void sflow_breaking_kb::tracking(lexer *p, fdm2D *b, ghostcell *pgc, slice &eta)
{
    if(p->count<=bt_count)
    return;

    bt_count = p->count;

    const double t = p->simtime;

    if(bt==nullptr)
    {
        bt = new breaking_crest3d(p);
        bt->Bon = p->A359;
        bt->slope_end = tan(p->A360*PI/180.0);
        bt->deep = p->A367;
        bt->Tdom = (p->B90>0) ? ((p->wTp>0.0) ? p->wTp : p->wT) : 0.0;
    }

    bt->Hmin = 2.0*p->A244;

    // the surface velocity: the depth-averaged velocity (uniform over the depth with the quadratic pressure)
    SLICELOOP4
    {
    us(i,j) = (p->wet[IJ]==1) ? b->U(i,j) : 0.0;
    vs(i,j) = (p->wet[IJ]==1 && p->j_dir==1) ? b->V(i,j) : 0.0;
    }

    pgc->gcsl_start4(p,us,1);
    pgc->gcsl_start4(p,vs,1);

    bt->detect(p,pgc,eta,us,vs,b->WL,zones_built ? brkzone.V : nullptr,t);

    // the seed: front faces and the band of A 502 L_k behind the breaking crest points, with their wavenumber
    auto sback = [&](const bc3_point &R) { return p->A502*((R.kw>0.0) ? tanh(R.kw*R.wl)/R.kw : R.wl); };
    auto integrand = [&](const bc3_point &R, double x, double y) { return 1.0; };
    auto epsf = [&](const bc3_point &R) { return 1.0; };

    SLICELOOP4
    {
    vbs(i,j) = 0.0;
    brkkseed(i,j) = 0.0;
    brkseed(i,j) = 0;
    }

    bt->rays(p,pgc,eta,b->WL,t,sback,integrand,epsf,vbs,brkseed.V,&brkkseed);

    pgc->gcsl_start4int(p,brkseed,50);
    pgc->gcsl_start4(p,brkkseed,1);

    if(p->P59==1)
    {
        if(p->mpirank==0)
        mkdir("./REEF3D_SFLOW_Breaking_Log",0777);

        bt->log(p,t,"./REEF3D_SFLOW_Breaking_Log/REEF3D-SFLOW-Breaking-Crests.dat");
    }

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking crest points: "<<bt->nbreak<<endl;
}

void sflow_breaking_kb::viscosity(lexer *p, fdm2D *b, ghostcell *pgc, slice &eta, slice &eta_n, double alpha)
{
    // nhflow_breaking::kennedy_part with A 550 4 (seeded), depth averaged
    const double g = fabs(p->W22);
    const double dts = MAX(alpha*p->dt,1.0e-20);

    // new time step: the ages and the wavenumbers of the events at the end of the last step
    if(p->count>count_n)
    {
    SLICELOOP4
    {
    brkage_n(i,j) = brkage(i,j);
    brkk_n(i,j) = brkk(i,j);
    brkr_n(i,j) = brkr(i,j);
    }

    pgc->gcsl_start4(p,brkage_n,1);
    pgc->gcsl_start4(p,brkk_n,1);
    pgc->gcsl_start4(p,brkr_n,1);

    count_n = p->count;
    }

    int count = 0;

    SLICELOOP4
    {
    brkage(i,j) = 0.0;
    brkk(i,j) = 0.0;
    brkr(i,j) = 0.0;
    brkflag(i,j) = 0;
    b->vb(i,j) = 0.0;

        if(p->wet[IJ]==1 && brkzone(i,j)==0 && p->count>1)
        {
        const double h = MAX(b->WL(i,j),0.0);
        const double et = (eta(i,j) - eta_n(i,j))/dts;

        // the largest age of the column and its neighbours, with the wavenumber and initial threshold of that event
        double age = brkage_n(i,j);
        double kw = brkk_n(i,j);
        double r0 = brkr_n(i,j);

        auto take = [&](int ii, int jj) { if(brkage_n(ii,jj)>age) { age = brkage_n(ii,jj); kw = brkk_n(ii,jj); r0 = brkr_n(ii,jj); } };
        take(i-1,j);
        take(i+1,j);

        if(p->j_dir==1)
        {
        take(i,j-1);
        take(i,j+1);
        }

        // a seeded column starts or continues an event with the wavenumber of the breaking crest; columns without
        // seed start events at deta/dt > A 554 c (as A 550 3)
        bool seedc = false;

        if(brkseed(i,j)==1)
        {
            if(age<=0.0 || kw<=0.0)
            kw = brkkseed(i,j);

            seedc = true;
        }

        const double Lk = (kw>0.0) ? tanh(kw*h)/kw : h;
        const double c = sqrt(g*Lk);
        const double tstar = p->A503*sqrt(Lk/g);

        if(age<=0.0 || r0<=0.0)
        r0 = (seedc && c>0.0) ? MAX(p->A554, et/c) : p->A554;

        const double etI = r0*c;
        const double etF = p->A501*c;
        const double ets = (age>0.0 && tstar>0.0) ? etI + MIN(age/tstar,1.0)*(etF - etI) : etI;

            if(et>ets && ets>0.0)
            {
            const double B = (et>=2.0*ets) ? 1.0 : et/ets - 1.0;

            b->vb(i,j) = B*p->A504*p->A504*Lk*et;
            brkage(i,j) = age + p->dt;
            brkk(i,j) = kw;
            brkr(i,j) = r0;
            brkflag(i,j) = 1;
            ++count;
            }

            // the event of a seeded column ages also below the threshold
            else
            if(seedc)
            {
            brkage(i,j) = age + p->dt;
            brkk(i,j) = kw;
            brkr(i,j) = r0;
            }
        }
    }

    // spread over the band of A 502 L_k (bd: the largest width in cells)
    int nband = 0;

    SLICELOOP4
    {
    double dxl = p->DXN[IP];

    if(p->j_dir==1)
    dxl = MIN(dxl,p->DYN[JP]);

    bd(i,j) = (p->wet[IJ]==1 && brkzone(i,j)==0) ? int(ceil(p->A502*MAX(b->WL(i,j),0.0)/dxl - 1.0e-9)) : 0;
    nband = MAX(nband,bd(i,j));

    vbk(i,j) = brkk(i,j);
    }

    nband = pgc->globalimax(nband);

    for(int qn=0; qn<nband; ++qn)
    {
    pgc->gcsl_start4(p,b->vb,1);
    pgc->gcsl_start4(p,vbk,1);

    SLICELOOP4
    {
    vbp(i,j) = b->vb(i,j);
    vbs(i,j) = vbk(i,j);

        if(qn<bd(i,j))
        {
        const double hb = MAX(b->WL(i,j),0.0);
        const double dxb = (p->j_dir==1) ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];
        auto inband = [&](double kv) { return !(kv>0.0) || qn < int(ceil(p->A502*tanh(kv*hb)/kv/dxb - 1.0e-9)); };
        auto take = [&](int ii, int jj) { if(b->vb(ii,jj)>vbp(i,j) && inband(vbk(ii,jj))) { vbp(i,j) = b->vb(ii,jj); vbs(i,j) = vbk(ii,jj); } };
        take(i-1,j);
        take(i+1,j);

        if(p->j_dir==1)
        {
        take(i,j-1);
        take(i,j+1);
        }
        }
    }

    SLICELOOP4
    {
    b->vb(i,j) = (p->wet[IJ]==1) ? vbp(i,j) : 0.0;
    vbk(i,j) = vbs(i,j);
    }
    }

    SLICELOOP4
    b->breaking_print(i,j) = double(brkflag(i,j));

    pgc->gcsl_start4(p,b->vb,1);
    pgc->gcsl_start4(p,b->breaking_print,1);

    count = pgc->globalisum(count);

    double vbmax = 0.0;
    SLICELOOP4
    vbmax = MAX(vbmax,b->vb(i,j));
    vbmax = pgc->globalmax(vbmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<"  vb max: "<<vbmax<<endl;
}
