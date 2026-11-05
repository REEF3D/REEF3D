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
#include"cfd_amr_cpm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4.h"
#include"turbulence.h"
#include"rans_io.h"
#include"sediment_part.h"
#include"reinitopo_void.h"
#include"reinitopo_RK3.h"
#include"ioflow_void.h"
#include"solver_void.h"
#include"patchBC_void.h"
#include"gridfile_v2.h"
#include"dmpart_api.h"
#include"grid_helper.h"
#include"mgcslice1.h"
#include"mgcslice2.h"
#include"mgcslice4.h"
#include<cmath>
#include<iostream>
#include<algorithm>

/*--------------------------------------------------------------------
CPM on the CFD mesh refinement, step 1 of the design note (validation 25): static refinement boxes
(G 15), one time step for all levels (G 7 0), the bed and all parcels within the finest level.

The grid of CPM. The cells of the finest level L whose centres lie in the bounding box of the boxes
G 15 make a regular grid. Its global grid file (the format of DIVEMesh M 50, gridfile_v2.h: the
nodes of level L, all cells fluid, the sides walls, periodic where level 0 is periodic and the box
spans the domain) is built in memory and decomposed over all ranks by dmpart, as REEF3D M 40 does
for level 0; every rank builds the lexer (lexer(const lexer&, int), the grid read from memory),
the ghostcell (its own communicator) and the fdm of its piece, and on them the objects of the
sediment step of a CFD run: sediment_part (CPM), the reinitialisation of the bed, an ioflow and a
solver without function, and a turbulence model that holds k and its dissipation (copied).

The copy. Every cell of the grid of CPM with its margin (global index of level L, wrapped at the
periodic sides) lies in a patch of level L (the boxes G 15 are refined to level L, cpm_setup
checks that the grid is covered); the items are requested once from the ranks of those patches
(MPI_Alltoallv), and every step the fluid goes over them onto the grid of CPM: u, v, w, the pressure,
the level set, the density, the viscosity, the eddy viscosity and k, its dissipation (RANS); at the
start also the solids and the bed. The bed (topo) of the cells inside the grid goes back onto the
patches; it is restricted onto the coarser levels (the mean of the children, cfd_amr::restrict_level)
down to level 0, the cells around the patches are filled from it, and every grid updates its direct
forcing as level 0 does after a bed change (the flags of D 22 1, gcdf_update).

The step. At the start of the time step, where level 0 runs its sediment step: the time and the
run state of level 0 onto the lexer of CPM (amr_follow), the fluid onto the grid, the sediment step
of CPM there (sediment_part::start_cfd), the bed back. The fluid step of the refinement follows.

Outside the grid of CPM the bed is fixed (as Q 76, Q 77 on one grid); a parcel at a side of the grid
that is not periodic meets a wall. The bed level of the columns of every grid (fdm::bed, read by the
wall laws of k-omega and k-epsilon with sediment) follows the bed. Not yet: the two-way coupling (Q 50 1),
the Eulerian suspension of the hybrid layer (Q 58 3, 4), subcycling (G 7 1), adaptive refinement, the
placement of the patches (G 40 1, 2), checkpoints and the hotstart (P 48, I 45, I 40), floating bodies.
The outputs of CPM (the parcels Q 180, the log of CPM) come from the grid of CPM; level 0 writes no bed
of its own.
--------------------------------------------------------------------*/

namespace
{
// a patch kernel: ghostcell without exchange, rank-local reductions, the fdm of the patch (as cfd_amr.cpp)
struct pscope
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};

// k and its dissipation of the patches, held for the bed shear stress of the layer (Q 68)
class turbulence_amr_copy : public turbulence
{
public:
    turbulence_amr_copy(lexer *p) : kin(p), eps(p) {}
    field4 kin, eps;

    void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*) override {}
    void ktimesave(lexer*, fdm*, ghostcell*) override {}
    void etimesave(lexer*, fdm*, ghostcell*) override {}
    void isource(lexer*, fdm*) override {}
    void jsource(lexer*, fdm*) override {}
    void ksource(lexer*, fdm*) override {}
    void print_3D(lexer*, fdm*, ghostcell*, std::vector<char>&, size_t&) override {}
    void ini(lexer*, fdm*, ghostcell*) override {}
    double kinval(int i, int j, int k) override { return kin(i,j,k); }
    double epsval(int i, int j, int k) override { return eps(i,j,k); }
    void gcupdate(lexer*, fdm*, ghostcell*) override {}
    double ccipol_kinval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4_kin(kin,x,y,z); }
    double ccipol_epsval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4(eps,x,y,z); }
    double ccipol_a_kinval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4_kin(kin,x,y,z); }
    double ccipol_a_epsval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4a(eps,x,y,z); }
    void kinget(int i, int j, int k, double v) override { kin(i,j,k) = v; }
    void epsget(int i, int j, int k, double v) override { eps(i,j,k) = v; }
    void name_ParaView_parallel(lexer*, ofstream&) override {}
    void name_ParaView(lexer*, std::stringstream&, int*, int&) override {}
    void offset_ParaView(lexer*, int*, int&) override {}
};
}

// --------------------------------------------------------------------- scope
// CPM (S 10 1, Q 10 >= 1) on the refinement: what this version covers; nullptr if it is covered
const char* cfd_amr::cpm_why(lexer *p)
{
    if(p->S10!=1 || p->Q10<1)
    return "no sediment, or CPM with the bed as topography (S 10 1, Q 10 1)";
    if(p->G15<1)
    return "CPM: the refinement boxes G 15 over the bed (static refinement)";
    if(p->G7!=0)
    return "CPM: one time step for all levels (G 7 0)";
    if(p->G40!=0)
    return "CPM: the patches on the ranks of the level-0 cells below (G 40 0)";
    if(p->G12>0 || p->G24>0 || p->G25>0 || p->G26>0)
    return "CPM: static refinement (no G 12, G 24, G 25, G 26)";
    if(p->Q50==1)
    return "CPM: no two-way coupling (Q 50 0)";
    if(p->Q58>=3)
    return "CPM: the bedload layer without the Eulerian suspension (Q 58 0-2)";
    if(p->P48>0 || p->P48_dt>0.0 || p->I45==1 || p->I40==1)
    return "CPM: no checkpoints or hotstart (P 48, I 45, I 40)";
    if(p->X10>0)
    return "CPM: no floating bodies (X 10 0)";
    return nullptr;
}

// --------------------------------------------------------------------- the grid of CPM
// the global grid file of the box (format of DIVEMesh M 50), rank 0
static std::vector<char> cpm_box_grid(lexer *p0, const int *n, const std::vector<double> *X, int px, int py, int nranks)
{
    const int xma=3, yma=3, zma=3, marge=increment::marge;
    const int NI=n[0]+2*xma, NJ=n[1]+2*yma, NK=n[2]+2*zma;
    auto c3 = [&](int i, int j, int k) { return ((size_t)(i+xma)*NJ + (size_t)(j+yma))*NK + (size_t)(k+zma); };

    // the boundary types of the sides (DIVEMesh C 11-16): walls, symmetry in y of a 2D grid; periodic
    // sides (C 21, C 22) have no boundary surfaces
    const int C11=21, C14=21, C15=21, C16=21;
    const int C12 = (n[1]==1) ? 3 : 21, C13 = (n[1]==1) ? 3 : 21;
    const int C21 = px ? 1 : 0, C22 = py ? 1 : 0, C23 = 0;

    std::vector<int> flag((size_t)NI*NJ*NK,-21), solid((size_t)NI*NJ*NK,-10), fsl((size_t)NI*NJ,-10);
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    {
        fsl[(size_t)(i+xma)*NJ + (j+yma)] = 1;
        for(int k=0; k<n[2]; ++k)
        {
            flag[c3(i,j,k)] = 1;
            solid[c3(i,j,k)] = 1;
        }
    }
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    for(int k=0; k<n[2]; ++k)
    {
        if(i==0) flag[c3(i-1,j,k)] = -C11;
        if(i==n[0]-1) flag[c3(i+1,j,k)] = -C14;
        if(j==0) flag[c3(i,j-1,k)] = -C13;
        if(j==n[1]-1) flag[c3(i,j+1,k)] = -C12;
        if(k==0) flag[c3(i,j,k-1)] = -C15;
        if(k==n[2]-1) flag[c3(i,j,k+1)] = -C16;
    }

    // the boundary surfaces (DIVEMesh surface::makesurf)
    std::vector<int> surf;
    auto add = [&](int i, int j, int k, int side, int nb)
    {
        surf.push_back(i); surf.push_back(j); surf.push_back(k); surf.push_back(side); surf.push_back(std::abs(nb));
    };
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    for(int k=0; k<n[2]; ++k)
    {
        if(flag[c3(i-1,j,k)]<0 && (C21<=1 || i>0)) add(i,j,k,1,flag[c3(i-1,j,k)]);
        if(flag[c3(i+1,j,k)]<0 && (C21<=1 || i<n[0]-1)) add(i,j,k,4,flag[c3(i+1,j,k)]);
        if(flag[c3(i,j-1,k)]<0 && (C22<=1 || j>0)) add(i,j,k,3,flag[c3(i,j-1,k)]);
        if(flag[c3(i,j+1,k)]<0 && (C22<=1 || j<n[1]-1)) add(i,j,k,2,flag[c3(i,j+1,k)]);
        if(flag[c3(i,j,k-1)]<0 && (C23<=1 || k>0)) add(i,j,k,5,flag[c3(i,j,k-1)]);
        if(flag[c3(i,j,k+1)]<0 && (C23<=1 || k<n[2]-1)) add(i,j,k,6,flag[c3(i,j,k+1)]);
    }
    const int nsurf = int(surf.size()/5);

    gridv2::writer w;
    w.start(gridv2::magic_glob);

    size_t pos = w.begin("GHED");
    {
        // the decomposition: dmpart chooses (M 20 1), the cuts in x (and y in 3D)
        const int iv[] = {n[0], n[1], n[2], xma, yma, zma, marge,
                          C11, C12, C13, C14, C15, C16, C21, C22, C23,
                          p0->i_dir, p0->j_dir, p0->k_dir, 0, 0, 0,
                          1, (n[1]>1 ? 1 : 0), 0, 1, 0, 0, 0, 0, 0, nranks,
                          nsurf};
        const int nint = int(sizeof(iv)/sizeof(int));
        w.put_int(nint);
        for(int q=0; q<nint; ++q)
        w.put_int(iv[q]);

        const double dv[] = {X[0][marge+1]-X[0][marge], 0.0, 0.0, 0.0,
                             X[0][marge], X[1][marge], X[2][marge], X[0][marge+n[0]], X[1][marge+n[1]], X[2][marge+n[2]],
                             p0->global_orig_x, p0->global_orig_y, p0->alpha_grid, 0.0, 0.0, 0.0};
        const int ndbl = int(sizeof(dv)/sizeof(double));
        w.put_int(ndbl);
        for(int q=0; q<ndbl; ++q)
        w.put_double(dv[q]);
    }
    w.end(pos);

    pos = w.begin("GFLG");
    w.put_rle(flag.size(),[&](size_t q){ return flag[q]; });
    w.end(pos);

    pos = w.begin("GSOL");
    w.put_rle(solid.size(),[&](size_t q){ return solid[q]; });
    w.end(pos);

    pos = w.begin("GSLF");
    w.put_rle(fsl.size(),[&](size_t q){ return fsl[q]; });
    w.end(pos);

    pos = w.begin("NODE");
    for(int d=0; d<3; ++d)
    for(int q=0; q<n[d]+1+2*marge; ++q)
    w.put_double(X[d][q]);
    w.end(pos);

    pos = w.begin("SURF");
    w.put_table((long long)nsurf,5,[&](long long r, int c){ return surf[(size_t)r*5+c]; });
    w.end(pos);

    w.finish();
    return std::move(w.buf);
}

void cfd_amr::cpm_setup(lexer *p, fdm *a, ghostcell *pgc)
{
    cfd_amr_cpm *C = cpm = new cfd_amr_cpm;
    const int marge = increment::marge;
    C->L = maxlev;

    // the box: the bounding box of the boxes G 15, the cells of level L whose centre lies in it
    double bx0[3], bx1[3];
    bx0[0] = p->G15_xs[0]; bx1[0] = p->G15_xe[0];
    bx0[1] = p->G15_ys[0]; bx1[1] = p->G15_ye[0];
    bx0[2] = p->G15_zs[0]; bx1[2] = p->G15_ze[0];
    for(int b=1; b<p->G15; ++b)
    {
        bx0[0] = MIN(bx0[0],p->G15_xs[b]); bx1[0] = MAX(bx1[0],p->G15_xe[b]);
        bx0[1] = MIN(bx0[1],p->G15_ys[b]); bx1[1] = MAX(bx1[1],p->G15_ye[b]);
        bx0[2] = MIN(bx0[2],p->G15_zs[b]); bx1[2] = MAX(bx1[2],p->G15_ze[b]);
    }
    for(int d=0; d<3; ++d)
    {
        const int N = gn(C->L,d);
        int lo = N, hi = -1;
        if(rr[d]==1)
        {
            lo = 0; hi = N-1;
        }
        else
        for(int I=0; I<N; ++I)
        {
            const double xc = 0.5*(node(C->L,d,I)+node(C->L,d,I+1));
            if(xc>=bx0[d] && xc<=bx1[d])
            {
                lo = MIN(lo,I);
                hi = MAX(hi,I);
            }
        }
        if(hi<lo)
        {
            if(myrank==0)
            cout<<endl<<"!!! CPM on the mesh refinement: the boxes G 15 hold no cell !!!"<<endl<<endl;
            pgc->final(true);
        }
        C->lo[d] = lo;
        C->n[d] = hi-lo+1;
    }

    // the grid covered by level L
    {
        long miss = 0;
        int I[3];
        for(I[0]=C->lo[0]; I[0]<C->lo[0]+C->n[0]; ++I[0])
        for(I[1]=C->lo[1]; I[1]<C->lo[1]+C->n[1]; ++I[1])
        for(I[2]=C->lo[2]; I[2]<C->lo[2]+C->n[2]; ++I[2])
        if(gpatch_at(C->L,I)<0)
        ++miss;
        if(miss>0)
        {
            if(myrank==0)
            cout<<endl<<"!!! CPM on the mesh refinement: "<<miss<<" cells of the box of the boxes G 15 are not on level "<<C->L<<" !!!"<<endl<<endl;
            pgc->final(true);
        }
    }

    // periodic where level 0 is periodic and the box spans the domain
    const int px = (per[0] && C->lo[0]==0 && C->n[0]==gn(C->L,0)) ? 1 : 0;
    const int py = (per[1] && C->lo[1]==0 && C->n[1]==gn(C->L,1)) ? 1 : 0;

    // the global grid file, decomposed by rank 0, the grid file of every rank sent to it
    vector<char> mine;
    {
        int ok = 1;
        if(myrank==0)
        {
            vector<double> X[3];
            for(int d=0; d<3; ++d)
            {
                X[d].resize(C->n[d]+1+2*marge);
                for(int q=0; q<(int)X[d].size(); ++q)
                X[d][q] = node(C->L,d,C->lo[d]+q-marge);
            }
            vector<char> glob = cpm_box_grid(p,C->n,X,px,py,nranks);
            vector<vector<char>> rb;
            string msg;
            if(!dmpart::build_mem(std::move(glob),nranks,rb,msg))
            {
                cout<<endl<<"!!! CPM on the mesh refinement: the decomposition of its grid: "<<msg<<" !!!"<<endl<<endl;
                ok = 0;
            }
            MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);
            if(ok)
            {
                for(int r=1; r<nranks; ++r)
                {
                    long long nb = (long long)rb[r].size();
                    MPI_Send(&nb,1,MPI_LONG_LONG,r,4731,MPI_COMM_WORLD);
                    MPI_Send(rb[r].data(),int(nb),MPI_CHAR,r,4732,MPI_COMM_WORLD);
                }
                mine = std::move(rb[0]);
            }
        }
        else
        {
            MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);
            if(ok)
            {
                long long nb = 0;
                MPI_Recv(&nb,1,MPI_LONG_LONG,0,4731,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                mine.resize(size_t(nb));
                MPI_Recv(mine.data(),int(nb),MPI_CHAR,0,4732,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        if(!ok)
        pgc->final(true);
    }

    // the lexer, ghostcell and fdm of the grid of CPM (the steps of driver for a grid of CFD)
    lexer *q = C->p = new lexer(*p,0);
    ghostcell *g = C->pgc = new ghostcell(q);
    q->grid_inmem = true;
    q->grid_mem = std::move(mine);
    q->read_grid();
    q->vellast();
    q->lexer_ini();
    q->amr_inherit(*p);
    q->solidread = p->solidread;
    q->toporead = p->toporead;
    q->topoforcing = p->topoforcing;
    q->G501 = 1;

    g->gc_ini(q);
    g->gcx_ini(q);
    q->gridini(g);

    q->flagini();
    g->flagfield(q);
    g->flagx(q,q->flagsf1);
    g->flagx(q,q->flagsf2);
    g->flagx(q,q->flagsf3);
    g->flagx(q,q->flagsf4);
    g->flagx(q,q->flag1);
    g->flagx(q,q->flag2);
    g->flagx(q,q->flag3);
    g->flagx(q,q->flag4);
    g->gcxupdate(q);
    q->vecsize(g);
    {
        grid_helper gh(q);
        gh.fillgcb1(q);
        gh.fillgcb2(q);
        gh.fillgcb3(q);
        gh.fillgcb4a(q);
        gh.fillgcb4_wall(q);
        gh.make_dgc(q);
        gh.fill_dgc1(q);
        gh.fill_dgc2(q);
        gh.fill_dgc3(q);
        gh.fill_dgc4(q);
    }
    {
        g->gcslflagx(q,q->flagslice4);
        mgcslice1 m1(q);
        mgcslice2 m2(q);
        mgcslice4 m4(q);
        m1.makemgc(q);
        g->gcslflagx(q,q->flagslice1);
        m1.gcb_seed(q);
        m2.makemgc(q);
        g->gcslflagx(q,q->flagslice2);
        m2.gcb_seed(q);
        m4.makemgc(q);
        m4.gcb_seed(q);
        g->gcsl_setbc1(q);
        g->gcsl_setbc2(q);
        g->gcsl_setbc4(q);
        g->gcsl_setbcio(q);
        g->dgcslini1(q);
        g->dgcslini2(q);
        g->dgcslini4(q);
    }

    C->a = new fdm(q);
    g->fdm_update(C->a);

    if(myrank==0)
    cout<<"CPM on the mesh refinement: the grid of level "<<C->L<<", "<<C->n[0]<<" x "<<C->n[1]<<" x "<<C->n[2]<<" cells from ("
        <<node(C->L,0,C->lo[0])<<", "<<node(C->L,1,C->lo[1])<<", "<<node(C->L,2,C->lo[2])<<") to ("
        <<node(C->L,0,C->lo[0]+C->n[0])<<", "<<node(C->L,1,C->lo[1]+C->n[1])<<", "<<node(C->L,2,C->lo[2]+C->n[2])<<")"
        <<(px ? ", periodic in x" : "")<<(py ? ", periodic in y" : "")<<endl;

    cpm_plans();

    // the objects of the sediment step
    C->pturb = new turbulence_amr_copy(q);
    C->pBC = new patchBC_void(q);
    C->pflow = new ioflow_v(q,g,C->pBC);
    C->psolv = new solver_void(q,C->a,g);
    if(p->G540>0)
    C->preto = new reinitopo_RK3(q);
    else
    C->preto = new reinitopo_void();

    // the fluid, the solids and the bed of the patches
    cpm_pull(true);

    g->solid_forcing_ini(q,C->a);

    C->psed = new sediment_part(q,C->a,g,C->pturb,C->pBC);
    C->psed->ini_cfd(q,C->a,g);
    C->psed->update_cfd(q,C->a,g,C->pflow,C->preto);

    cpm_push();
}

// --------------------------------------------------------------------- the copy plans
void cfd_amr::cpm_plans()
{
    cfd_amr_cpm *C = cpm;
    lexer *q = C->p;
    const int m = q->margin;      // the margin of the fields (3), not the nodes (marge)

    // the cells of this piece with their margin: the patch of level L that holds each
    vector<vector<int>> req(nranks);
    vector<vector<int>> pos(nranks);
    vector<vector<char>> inn(nranks);
    int I[3], J[3];
    for(int i=-m; i<q->knox+m; ++i)
    for(int j=(q->j_dir==1 ? -m : 0); j<(q->j_dir==1 ? q->knoy+m : 1); ++j)
    for(int k=-m; k<q->knoz+m; ++k)
    {
        I[0] = C->lo[0] + q->origin_i + i;
        I[1] = C->lo[1] + q->origin_j + j;
        I[2] = C->lo[2] + q->origin_k + k;
        if(!wrap(C->L,I,J))
        continue;
        const int gp = gpatch_at(C->L,J);
        if(gp<0)
        continue;
        const r3gpatch &G = GP[gp];
        req[G.rank].push_back(G.lid);
        req[G.rank].push_back(J[0]-G.lo[0]);
        req[G.rank].push_back(J[1]-G.lo[1]);
        req[G.rank].push_back(J[2]-G.lo[2]);
        req[G.rank].push_back((i>=0 && i<q->knox && j>=0 && j<q->knoy && k>=0 && k<q->knoz) ? 1 : 0);
        pos[G.rank].push_back((i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin);
        inn[G.rank].push_back((i>=0 && i<q->knox && j>=0 && j<q->knoy && k>=0 && k<q->knoz) ? 1 : 0);
    }

    C->scnt.assign(nranks,0);
    C->rcnt.assign(nranks,0);
    C->sdsp.assign(nranks,0);
    C->rdsp.assign(nranks,0);
    for(int r=0; r<nranks; ++r)
    C->scnt[r] = (int)pos[r].size();
    MPI_Alltoall(C->scnt.data(),1,MPI_INT,C->rcnt.data(),1,MPI_INT,MPI_COMM_WORLD);
    for(int r=1; r<nranks; ++r)
    {
        C->sdsp[r] = C->sdsp[r-1]+C->scnt[r-1];
        C->rdsp[r] = C->rdsp[r-1]+C->rcnt[r-1];
    }
    const int ns = C->sdsp[nranks-1]+C->scnt[nranks-1];
    const int nr = C->rdsp[nranks-1]+C->rcnt[nranks-1];

    C->mine.clear();
    C->inner.clear();
    vector<int> sb;
    for(int r=0; r<nranks; ++r)
    {
        C->mine.insert(C->mine.end(),pos[r].begin(),pos[r].end());
        C->inner.insert(C->inner.end(),inn[r].begin(),inn[r].end());
        sb.insert(sb.end(),req[r].begin(),req[r].end());
    }

    vector<int> sc5(nranks), sd5(nranks), rc5(nranks), rd5(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc5[r] = 5*C->scnt[r]; sd5[r] = 5*C->sdsp[r];
        rc5[r] = 5*C->rcnt[r]; rd5[r] = 5*C->rdsp[r];
    }
    vector<int> rb(5*(size_t)nr);
    MPI_Alltoallv(sb.data(),sc5.data(),sd5.data(),MPI_INT,rb.data(),rc5.data(),rd5.data(),MPI_INT,MPI_COMM_WORLD);

    // served: the patch and the flat index of the cell in its arrays
    C->sv_g.resize(nr);
    C->sv_q.resize(nr);
    C->sv_in.resize(nr);
    for(int t=0; t<nr; ++t)
    {
        const int id = rb[5*t];
        lexer *pp = P[id]->pp;
        const int i = rb[5*t+1], j = rb[5*t+2], k = rb[5*t+3];
        C->sv_g[t] = id;
        C->sv_q[t] = (i-pp->imin)*pp->jmax*pp->kmax + (j-pp->jmin)*pp->kmax + k-pp->kmin;
        C->sv_in[t] = (char)rb[5*t+4];
    }
    (void)ns;
}

// the fluid of the patches onto the grid of CPM (ini: also the solids and the bed)
void cfd_amr::cpm_pull(bool ini)
{
    cfd_amr_cpm *C = cpm;
    fdm *a = C->a;
    turbulence_amr_copy *T = static_cast<turbulence_amr_copy*>(C->pturb);

    vector<std::function<field&(int)>> src;
    vector<field*> dst;
    auto add = [&](std::function<field&(int)> s, field &d) { src.push_back(s); dst.push_back(&d); };
    add([this](int g) -> field& { return gfd(g)->u; }, a->u);
    add([this](int g) -> field& { return gfd(g)->v; }, a->v);
    add([this](int g) -> field& { return gfd(g)->w; }, a->w);
    add([this](int g) -> field& { return gfd(g)->press; }, a->press);
    add([this](int g) -> field& { return gfd(g)->phi; }, a->phi);
    add([this](int g) -> field& { return gfd(g)->ro; }, a->ro);
    add([this](int g) -> field& { return gfd(g)->visc; }, a->visc);
    add([this](int g) -> field& { return gfd(g)->eddyv; }, a->eddyv);
    if(turb)
    {
        add([this](int g) -> field& { return kfield(g,0); }, T->kin);
        add([this](int g) -> field& { return kfield(g,1); }, T->eps);
    }
    if(ini)
    {
        add([this](int g) -> field& { return gfd(g)->topo; }, a->topo);
        add([this](int g) -> field& { return gfd(g)->solid; }, a->solid);
    }
    const int nf = (int)src.size();

    // the served items, field by field per rank (the order of the requests)
    const int nr = (int)C->sv_g.size();
    vector<double> sb((size_t)nr*nf);
    for(int r=0, t=0; r<nranks; ++r)
    for(int u=0; u<C->rcnt[r]; ++u, ++t)
    for(int f=0; f<nf; ++f)
    sb[(size_t)t*nf+f] = src[f](C->sv_g[t]).V[C->sv_q[t]];

    vector<int> sc(nranks), sd(nranks), rc(nranks), rd(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc[r] = nf*C->rcnt[r]; sd[r] = nf*C->rdsp[r];
        rc[r] = nf*C->scnt[r]; rd[r] = nf*C->sdsp[r];
    }
    const size_t nm = C->mine.size();
    vector<double> rb(nm*nf);
    MPI_Alltoallv(sb.data(),sc.data(),sd.data(),MPI_DOUBLE,rb.data(),rc.data(),rd.data(),MPI_DOUBLE,MPI_COMM_WORLD);

    for(size_t t=0; t<nm; ++t)
    for(int f=0; f<nf; ++f)
    dst[f]->V[C->mine[t]] = rb[t*nf+f];
}

// the bed of the grid of CPM back onto the patches of level L, restricted to level 0, the cells around
// the patches filled, the direct forcing of all grids
void cfd_amr::cpm_push()
{
    cfd_amr_cpm *C = cpm;
    fdm *a = C->a;

    const size_t nm = C->mine.size();
    vector<double> sb(nm);
    for(size_t t=0; t<nm; ++t)
    sb[t] = a->topo.V[C->mine[t]];

    const int nr = (int)C->sv_g.size();
    vector<double> rb(nr);
    MPI_Alltoallv(sb.data(),C->scnt.data(),C->sdsp.data(),MPI_DOUBLE,rb.data(),C->rcnt.data(),C->rdsp.data(),MPI_DOUBLE,MPI_COMM_WORLD);

    for(int t=0; t<nr; ++t)
    if(C->sv_in[t])
    CP(C->sv_g[t])->a->topo.V[C->sv_q[t]] = rb[t];

    // restriction: level L down to level 0
    vector<std::function<field&(int)>> cells = {[this](int g) -> field& { return gfd(g)->topo; }};
    for(int l=C->L; l>=1; --l)
    restrict_level(l,cells,{});

    // level 0: its halo, the direct forcing
    pgc0->start4a(p0,a0->topo,150);
    if(p0->D22==1)
    pgc0->solid_forcing_flag_update(p0,a0);
    pgc0->gcdf_update(p0,a0);

    // the cells around the patches, level by level; the direct forcing of the patches
    vector<fspec> fs(1);
    fs[0].type = CELL;
    fs[0].mode = 0;
    fs[0].f = [this](int g) -> field& { return gfd(g)->topo; };
    for(int l=1; l<=maxlev; ++l)
    fill(l,fs);

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope ps(pgc0,c->a,a0);
        if(p0->D22==1)
        pgc0->solid_forcing_flag_update(c->pp,c->a);
        pgc0->gcdf_update(c->pp,c->a);
    }

    // the bed level of the columns (the wall law of k-omega and k-epsilon samples the velocity above it
    // with sediment, S 10 > 0): from the bed of every grid
    cpm_bed(p0,a0);
    pgc0->gcsl_start4(p0,a0->bed,50);
    for(int id=0; id<npat; ++id)
    cpm_bed(CP(id)->pp,CP(id)->a);
}

// the bed level of the columns of a grid from its bed (topo): where it crosses zero (as CPM::bedzh_update),
// the bottom of the domain where the column has no bed; a column of a patch whose bed lies below it keeps
// its level
void cfd_amr::cpm_bed(lexer *pp, fdm *fa)
{
    const int m = increment::marge;
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    {
        bool set = false;
        double h = pp->ZN[0+m];
        if(fa->topo(i,j,0)>=0.0 && pp->origin_k==0)
        set = true;
        for(int k=1; k<pp->knoz; ++k)
        if(fa->topo(i,j,k-1)<0.0 && fa->topo(i,j,k)>=0.0)
        {
            h = -(fa->topo(i,j,k-1)*pp->DZP[k-1+m])/(fa->topo(i,j,k)-fa->topo(i,j,k-1)) + pp->ZP[k+m]-pp->DZP[k-1+m];
            set = true;
        }
        if(set)
        fa->bed(i,j) = h;
    }
}

// --------------------------------------------------------------------- the step
void cfd_amr::cpm_step(lexer *p, fdm *a, ghostcell *pgc)
{
    if(cpm==nullptr)
    return;

    cfd_amr_cpm *C = cpm;
    lexer *q = C->p;

    // the time of level 0, and the wave period (the start of the sediment S 41 3); the sediment times
    // (sedtime, ...) are the grid's own here, not the run state of level 0 (amr_follow)
    q->wT = p->wT;
    q->simtime = p->simtime;
    q->dt = p->dt;
    q->dt_old = p->dt_old;
    q->count = p->count;

    cpm_pull(false);
    C->psed->start_cfd(q,C->a,C->pgc,C->pflow,C->preto,C->psolv);
    cpm_push();
}
