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

#include"nhflow_amr.h"
#include"nhflow_amr_fill.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"nhflow_momentum_func.h"
#include"nhflow_fsf.h"
#include"nhflow_potential_f.h"
#include"ioflow_amr_block.h"
#include<cstdlib>
#include<cstring>
#include<cstdint>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>

// Multi-block level 0 of REEF3D::NHFLOW (G 50 1 with the blocks of DIVEMesh R 20 1, phase P4 of the
// multi-block curvilinear concept), the glue of REEFAMR (reefamr_block.cpp) as for SFLOW:
//  - block 0 is the native grid: its own cells and a ring of cells on every face with an interface
//    (the cells of the neighbour blocks); blocks 1.. are level-0 patches (lev[0]) with their own
//    lexer (the block nodes), fdm_nhf and NHFLOW objects, placed on the ranks for the load (G 55: the
//    blocks of DIVEMesh cut into pieces first, every piece a block, geo_mesh_split.cpp)
//  - before every stage the cells around the blocks and the ring take the stage input of the cells
//    they stand for (eta, water level, wet and deep, U, V, W, UH, VH, WH per layer, P per node):
//    the stage of every block then runs as on one grid, velocities are Cartesian components, so the
//    values need no rotation; the bed of these cells is the bed of their source (fill_bed0)
//  - the non-hydrostatic pressure is one composite solve over the cells of all blocks (and the
//    patches): the ring and the cells around the blocks are filled from the vectors, the
//    preconditioner runs a V-cycle per block (nhflow_amr_press.cpp)
//  - faces of the blocks 1.. without an interface are walls, except the open faces of G 54 (discharge
//    or water level, constant or a time series) and the zones of B 522 (Riemann edges, with the waves of
//    B 524 sources, Flather, clamped level or discharge edges with a background; relaxation and beach
//    zones) on their -i and +i faces (ioflow_amr_block; a block with such faces on its j faces is turned,
//    geo_mesh::turn_blocks)
//  - breaking (A 550 1): block 0 and the blocks run the parts of nhflow_breaking together, the ring and the
//    cells around the blocks take the detection, the flags and the viscosity of the cells they stand for
//    (breaking_level0, nhflow_amr_diff.cpp)

using namespace nhflow_amr_detail;

namespace
{
typedef reefamr_comms_off comms_off;
}

// what runs on a multi-block grid
void nhflow_amr::blocks_check(lexer *p)
{
    int ok=1;
    if(p->G7==1 && p->X10>0) ok=0;  // subcycling (sub_level0_blocks): without bodies
    if(p->X10>2) ok=0;              // floating (X 10 1) and prescribed (X 10 2) bodies only (P5c)
    if(p->X10>0 && p->G40>0) ok=0;  // bodies: the patches on the rank of the block piece below
    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW multi-block grid: not with subcycling (G 7 1) and bodies, the ship pressure X 10 3, or X 10 1/2 with G 40 !!!"<<endl;
        exit(1);
    }
}

// blocks 1..: bed from the block nodes (mean of the four node levels, as lexer::geometry_bed), state
// as nhflow_f::ini and nhflow_fsf_f::ini on level 0 (eta 0, the F 72 boxes by the cell centre, at rest)
void nhflow_amr::blocks_ini(lexer *p, ghostcell *pgc)
{
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        fdm_nhf *d = c->d;
        pscope ps(pgc,d,d0);
        c->fresh = false;

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            const int i = ii-EXT, j = jj-EXT;
            double x,y,z00,z10,z01,z11;
            block_node_any(c->blk,i,j,x,y,z00);
            block_node_any(c->blk,i+1,j,x,y,z10);
            block_node_any(c->blk,i,j+1,x,y,z01);
            block_node_any(c->blk,i+1,j+1,x,y,z11);
            const double zb = 0.25*(z00 + z10 + z01 + z11);
            d->bed(ii,jj) = MAX(p->global_zmin,zb);
        }

        // the wall ghost cells as level 0 has them (gcsl_start4), then the smoothing sweeps of nhflow_f::ini (A 509), over the cells of the block and around it (the
        // cells around it then take the bed of their source, fill_bed0)
        // (G 52 1 on multi-block grids: from the old values, as nhflow_f::ini)
        pgc->gcsl_start4(pp,d->bed,50);
        for(int qn=0; qn<p->A509; ++qn)
        {
            vector<double> b0((size_t)pp->imax*pp->jmax);
            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            b0[lij(pp,ii,jj)] = d->bed(ii,jj);
            auto B = [&](int a, int b) { return b0[lij(pp,a,b)]; };
            for(int ii=0; ii<pp->knox; ++ii)
            for(int jj=0; jj<pp->knoy; ++jj)
            if(pp->flagslice4[lij(pp,ii,jj)]>0)
            d->bed(ii,jj) = 0.5*B(ii,jj) + 0.125*(B(ii-1,jj) + B(ii+1,jj) + B(ii,jj-1) + B(ii,jj+1));
            pgc->gcsl_start4(pp,d->bed,50);
        }

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            d->depth(ii,jj) = pp->wd - d->bed(ii,jj);
            pp->bed[lij(pp,ii,jj)] = d->bed(ii,jj);
        }

        const int K = pp->knoz;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            const int q = lij(pp,ii,jj);
            double eta = 0.0;
            for(int qn=0; qn<p->F72; ++qn)
            if(pp->XC2D[q]>=p->F72_xs[qn] && pp->XC2D[q]<p->F72_xe[qn] && pp->YC2D[q]>=p->F72_ys[qn] && pp->YC2D[q]<p->F72_ye[qn])
            eta = p->F72_h[qn] - p->F60;

            d->eta(ii,jj) = eta;
            d->eta_n(ii,jj) = eta;
            d->detadt(ii,jj) = 0.0;
            if(pp->flagslice4[q]>0)
            {
            pp->wet[q] = pp->wet_n[q] = (eta + d->depth(ii,jj)<p->A544) ? 0 : 1;
            pp->deep[q] = 1;
            }
            d->WL(ii,jj) = MAX(p->A544,eta + d->depth(ii,jj));
            pp->WL[q] = d->WL(ii,jj);

            for(int kk=0; kk<K; ++kk)
            {
                const int n = cidx(pp,ii,jj,kk);
                d->U[n] = d->V[n] = d->W[n] = 0.0;
                d->UH[n] = d->VH[n] = d->WH[n] = 0.0;
            }
            for(int kk=0; kk<=K; ++kk)
            d->P[fidx(pp,ii,jj,kk)] = 0.0;
        }
    }
}

// sigma grid and derived fields of the blocks 1.. (after their beds and the fills), as regrid_state
void nhflow_amr::blocks_sigma(ghostcell *pgc)
{
    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        fdm_nhf *d = c->d;
        pscope ps(pgc,d,d0);

        c->pmom->sigma_ini(pp,d,pgc,d->eta);
        c->pmom->inidisc(pp,d,pgc,c->pfsf);

        // the relaxation zones of the block (B 522) at the start (ioflow_amr_block::zones_ini), the sigma grid of the
        // surface they set
        if(c->pio!=nullptr && c->pio->has_zones())
        {
            c->pio->zones_ini(pp,d,pgc);
            c->pmom->sigma_ini(pp,d,pgc,d->eta);
            c->pmom->inidisc(pp,d,pgc,c->pfsf);
        }
        c->pfsf->kinematic_fsf(pp,d,d->U,d->V,d->W,d->eta);

        // floating bodies: the hull on the sigma grid of the block
        body_patch(pgc,*c);
    }
}

// the bed of the cells around the blocks and of the ring: the bed of the cell they stand for
void nhflow_amr::fill_bed0(ghostcell *pgc)
{
    fill_run0(1,7097,
              [&](const reefamr_fill &f, double *v)
              {
                  v[0] = (f.kind==0) ? gfd(f.g)->bed(f.si,f.sj) : 0.0;
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  lexer *pp = q->pp;
                  fdm_nhf *d = NP(q)->d;
                  d->bed(f.di,f.dj) = w[0];
                  d->depth(f.di,f.dj) = pp->wd - w[0];
                  pp->bed[lij(pp,f.di,f.dj)] = w[0];
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  d0->bed(f.di,f.dj) = w[0];
                  d0->depth(f.di,f.dj) = p0->wd - w[0];
                  p0->bed[lij(p0,f.di,f.dj)] = w[0];
              });

    {
    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        pgc->gcsl_start4(c->pp,c->d->bed,50);
        pgc->gcsl_start4(c->pp,c->d->depth,50);
    }
    }

    pgc->gcsl_start4(p0,d0->bed,50);
    pgc->gcsl_start4(p0,d0->depth,50);
}

// the cells around the blocks and the ring of the native grid with the input of stage s (s<0: the
// end of the step), copied from the cells they stand for
void nhflow_amr::fill_level0(ghostcell *pgc, int s)
{
    const int K = p0->knoz;
    const int nv = 4 + 6*K + K+1;

    fill_run0(nv,7098,
              [&](const reefamr_fill &f, double *v)
              {
                  if(f.kind!=0)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = 0.0;
                      return;
                  }
                  lexer *q = glex(f.g);
                  fdm_nhf *d = gfd(f.g);
                  stg S = stage_in(f.g,s);
                  v[0] = d->eta(f.si,f.sj);
                  v[1] = (*S.WL)(f.si,f.sj);
                  v[2] = q->wet[lij(q,f.si,f.sj)];
                  v[3] = q->deep[lij(q,f.si,f.sj)];
                  for(int k=0; k<K; ++k)
                  {
                      const int n = cidx(q,f.si,f.sj,k);
                      v[4+k] = d->U[n];
                      v[4+K+k] = d->V[n];
                      v[4+2*K+k] = d->W[n];
                      v[4+3*K+k] = S.UH[n];
                      v[4+4*K+k] = S.VH[n];
                      v[4+5*K+k] = S.WH[n];
                  }
                  for(int k=0; k<=K; ++k)
                  v[4+6*K+k] = d->P[fidx(q,f.si,f.sj,k)];
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  lexer *pp = q->pp;
                  fdm_nhf *d = NP(q)->d;
                  stg S = stage_in(id,s);
                  const int ii=f.di, jj=f.dj;
                  d->eta(ii,jj) = w[0];
                  (*S.WL)(ii,jj) = w[1];
                  pp->wet[lij(pp,ii,jj)] = (int)w[2];
                  pp->deep[lij(pp,ii,jj)] = (int)w[3];
                  for(int k=0; k<K; ++k)
                  {
                      const int n = cidx(pp,ii,jj,k);
                      d->U[n] = w[4+k];
                      d->V[n] = w[4+K+k];
                      d->W[n] = w[4+2*K+k];
                      S.UH[n] = w[4+3*K+k];
                      S.VH[n] = w[4+4*K+k];
                      S.WH[n] = w[4+5*K+k];
                  }
                  for(int k=0; k<=K; ++k)
                  d->P[fidx(pp,ii,jj,k)] = w[4+6*K+k];
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  lexer *q = p0;
                  fdm_nhf *d = d0;
                  stg S = stage_in(-1,s);
                  const int ii=f.di, jj=f.dj;
                  d->eta(ii,jj) = w[0];
                  (*S.WL)(ii,jj) = w[1];
                  q->wet[lij(q,ii,jj)] = (int)w[2];
                  q->deep[lij(q,ii,jj)] = (int)w[3];
                  for(int k=0; k<K; ++k)
                  {
                      const int n = cidx(q,ii,jj,k);
                      d->U[n] = w[4+k];
                      d->V[n] = w[4+K+k];
                      d->W[n] = w[4+2*K+k];
                      S.UH[n] = w[4+3*K+k];
                      S.VH[n] = w[4+4*K+k];
                      S.WH[n] = w[4+5*K+k];
                  }
                  for(int k=0; k<=K; ++k)
                  d->P[fidx(q,ii,jj,k)] = w[4+6*K+k];
              });

    {
    comms_off guard(pgc);
    for(int n : lev[0])
    bc_patch(pgc,*NP(n),s);
    }

    ring_halo(pgc,s);
}

// the ring cells in the partition halos of the native grid (exchange only: the boundary conditions
// of level 0 are left as they are)
void nhflow_amr::ring_halo(ghostcell *pgc, int s)
{
    if(p0->mpi_size==1)
    return;

    stg S = stage_in(-1,s);
    pgc->gcslparax(p0,d0->eta,4);
    pgc->gcslparacox(p0,d0->eta,gcval_eta);
    pgc->gcslparax(p0,*S.WL,4);
    pgc->gcslparacox(p0,*S.WL,gcval_eta);
    pgc->gcslparaxV_int(p0,p0->wet,4);
    pgc->gcslparacoxV_int(p0,p0->wet,50);
    pgc->gcslparaxV_int(p0,p0->deep,4);
    pgc->gcslparacoxV_int(p0,p0->deep,50);
    double *v[6] = {d0->U,d0->V,d0->W,S.UH,S.VH,S.WH};
    const int gv[6] = {10,11,12,14,15,16};
    for(int m=0; m<6; ++m)
    {
        pgc->gcparaxV(p0,v[m],gv[m]);
        pgc->gcparacoxV(p0,v[m],gv[m]);
    }
    pgc->gcparax7(p0,d0->P,7);
    pgc->gcparax7co(p0,d0->P,7);
}

// continuity (phase 0) or momentum (phase 1) part of stage s on the blocks 1..
void nhflow_amr::blocks_phase(ghostcell *pgc, int s, int phase)
{
    comms_off guard(pgc);
    for(int n : lev[0])
    {
        nhflow_amr_patch *c = NP(n);
        pscope ps(pgc,c->d,d0);
        if(phase==0)
        c->pmom->phase_F(c->pp,c->d,pgc,c->S,s);
        else
        c->pmom->phase_M(c->pp,c->d,pgc,c->S,s);
    }
}

// tests (environment REEF3D_AMR_DUMP=<dir>): the state of every grid of the rank (own cells), in the
// format of regression_dump, overwritten at every call
void nhflow_amr::dump_grids(lexer *p)
{
    static const char *dir = std::getenv("REEF3D_AMR_DUMP");
    if(dir==nullptr)
    return;

    vector<string> names;
    vector<vector<double>> data;

    auto grid = [&](int id, lexer *q, fdm_nhf *d, int lv, int b, int I0, int J0, int i0, int nx, int ny)
    {
        auto add = [&](const char *v)
        {
            std::ostringstream s;
            s<<"g"<<id<<"_"<<v;
            names.push_back(s.str());
            data.push_back(vector<double>());
            return &data.back();
        };
        vector<double> *a = add("box");
        a->push_back(lv); a->push_back(b); a->push_back(I0); a->push_back(J0); a->push_back(nx); a->push_back(ny); a->push_back(q->knoz);

        const char *nm2[3] = {"eta","WL","bed"};
        for(int m=0; m<3; ++m)
        {
            vector<double> *f = add(nm2[m]);
            for(int ii=i0; ii<i0+nx; ++ii)
            for(int jj=i0; jj<i0+ny; ++jj)
            f->push_back(m==0 ? d->eta(ii,jj) : (m==1 ? d->WL(ii,jj) : d->bed(ii,jj)));
        }
        const char *nm3[6] = {"U","V","W","UH","VH","WH"};
        double *src[6] = {d->U,d->V,d->W,d->UH,d->VH,d->WH};
        for(int m=0; m<6; ++m)
        {
            vector<double> *f = add(nm3[m]);
            for(int ii=i0; ii<i0+nx; ++ii)
            for(int jj=i0; jj<i0+ny; ++jj)
            for(int kk=0; kk<q->knoz; ++kk)
            f->push_back(src[m][cidx(q,ii,jj,kk)]);
        }
        vector<double> *f = add("P");
        for(int ii=i0; ii<i0+nx; ++ii)
        for(int jj=i0; jj<i0+ny; ++jj)
        for(int kk=0; kk<=q->knoz; ++kk)
        f->push_back(d->P[fidx(q,ii,jj,kk)]);

        // A 550: the breaking viscosity
        if(p->A550==1)
        {
            vector<double> *fv = add("vb");
            for(int ii=i0; ii<i0+nx; ++ii)
            for(int jj=i0; jj<i0+ny; ++jj)
            fv->push_back(d->vb(ii,jj));
        }

        // A 560: the eddy viscosity and k
        // (one field after the other: add() may move the vectors of the fields before)
        if(turb)
        {
            const char *nmt[2] = {"EV","KIN"};
            double *srt[2] = {d->EV,d->KIN};
            for(int m=0; m<2; ++m)
            {
                vector<double> *f = add(nmt[m]);
                for(int ii=i0; ii<i0+nx; ++ii)
                for(int jj=i0; jj<i0+ny; ++jj)
                for(int kk=0; kk<q->knoz; ++kk)
                f->push_back(srt[m][cidx(q,ii,jj,kk)]);
            }
        }

        // floating bodies: the level set of the hull and the solid flags
        if(b6!=nullptr)
        {
            vector<double> *fb = add("FB");
            for(int ii=i0; ii<i0+nx; ++ii)
            for(int jj=i0; jj<i0+ny; ++jj)
            for(int kk=0; kk<q->knoz; ++kk)
            fb->push_back(d->FB[cidx(q,ii,jj,kk)]);
            vector<double> *df = add("DF");
            for(int ii=i0; ii<i0+nx; ++ii)
            for(int jj=i0; jj<i0+ny; ++jj)
            for(int kk=0; kk<q->knoz; ++kk)
            df->push_back(q->DF[cidx(q,ii,jj,kk)]);
        }
    };

    grid(0,p0,d0,0,0,O0i,O0j,0,p0->knox,p0->knoy);
    for(int n=0; n<(int)P.size(); ++n)
    grid(n+1,P[n]->pp,NP(n)->d,P[n]->lev,P[n]->blk,P[n]->I0,P[n]->J0,EXT,P[n]->nx,P[n]->ny);

    std::ostringstream fn;
    fn<<dir<<"/amr_r"<<p->mpirank<<".bin";
    std::ofstream out(fn.str().c_str(),std::ios::binary);
    const char magic[8] = {'R','3','D','R','E','G','0','1'};
    out.write(magic,8);
    int32_t ival;
    double dval;
    ival=p->mpirank;  out.write((char*)&ival,sizeof(ival));
    ival=p->M10;      out.write((char*)&ival,sizeof(ival));
    ival=p->count;    out.write((char*)&ival,sizeof(ival));
    dval=p->simtime;  out.write((char*)&dval,sizeof(dval));
    ival=int32_t(names.size()); out.write((char*)&ival,sizeof(ival));
    for(size_t q=0; q<names.size(); ++q)
    {
        char name[16];
        std::memset(name,0,16);
        std::strncpy(name,names[q].c_str(),15);
        out.write(name,16);
        int64_t num = int64_t(data[q].size());
        out.write((char*)&num,sizeof(num));
        out.write((char*)data[q].data(),num*sizeof(double));
    }
}

// I 11 on a multi-block grid: the potential flow solve of nhflow_potential_f over the blocks, one
// composite solve as the pressure (the rows of every block, the ring and the cells around the blocks
// filled from the vectors, BiCGStab with the block-local V-cycles), then the velocities of every block
void nhflow_amr::potential(lexer *p, ghostcell *pgc)
{
    if(p->I11!=1)
    return;

    if(p->mpirank==0)
    cout<<"potential flow solver on "<<nblk<<" blocks..."<<endl;

    if(fabs(p->Ui)<1.0e-9 || fabs(p->Uo)<1.0e-9)
    {
        if(p->mpirank==0)
        cout<<"no inflow or outflow...potential flow stop"<<endl;
        return;
    }

    const double t0 = MPI_Wtime();
    nhflow_potential_f pot(p);
    auto len = [](lexer *q) { return (size_t)q->imax*q->jmax*(q->kmax+2); };

    vector<double> psi0(len(p0),0.0);
    vector<vector<double>> psi(P.size());

    pot.amr_assemble(p0,d0,pgc,psi0.data());
    ptgt0 = psi0.data();
    {
    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        psi[id].assign(len(c->pp),0.0);
        pot.amr_assemble(c->pp,c->d,pgc,psi[id].data());
        c->ptgt = psi[id].data();
    }
    }

    const int n46 = p->N46;
    p->N46 = 5000;
    pr_core(p,pgc);
    p->N46 = n46;
    const int it = pr_it_last;
    const double res = pr_res_last;
    --pr_solves;
    pr_it_total -= it;

    int fail = pot.amr_velocities(p0,d0,pgc,psi0.data());
    {
    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        fail = MAX(fail,pot.amr_velocities(c->pp,c->d,pgc,psi[id].data()));
    }
    }

    if(it>=5000 || !std::isfinite(res) || res>1.0e-2)
    fail = 1;
    fail = pgc->globalimax(fail);

    // as nhflow_potential_f::start: no converged solution (wet domain disconnected, in- and outflow
    // unbalanced): start from rest
    auto rest = [&](lexer *q, fdm_nhf *d)
    {
        const int n = q->imax*q->jmax*q->kmax;
        for(int m=0; m<n; ++m)
        d->U[m] = d->V[m] = d->W[m] = d->UH[m] = d->VH[m] = d->WH[m] = 0.0;
    };
    if(fail==1)
    {
        if(p->mpirank==0)
        cout<<"potential flow solver: no converged solution (wet domain disconnected or in/outflow unbalanced?) - starting from rest instead"<<endl;
        rest(p0,d0);
        for(int id : lev[0])
        rest(NP(id)->pp,NP(id)->d);
    }

    double *v[6] = {d0->U,d0->V,d0->W,d0->UH,d0->VH,d0->WH};
    const int gv[6] = {10,11,12,14,15,16};
    for(int m=0; m<6; ++m)
    pgc->start4V(p0,v[m],gv[m]);

    {
    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        double *w[6] = {c->d->U,c->d->V,c->d->W,c->d->UH,c->d->VH,c->d->WH};
        for(int m=0; m<6; ++m)
        pgc->start4V(c->pp,w[m],gv[m]);
    }
    }

    fill_level0(pgc,-1);

    if(p->mpirank==0)
    cout<<"potential flow on the blocks: iterations "<<it<<", residual "<<res<<", time "<<MPI_Wtime()-t0<<" s"<<endl<<endl;
}
