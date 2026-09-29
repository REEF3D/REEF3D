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

#include"sflow_amr.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice.h"
#include"sflow_HLL.h"
#include"sflow_signal_speed.h"
#include"sflow_reconstruct_hires.h"
#include"sflow_reconstruct_weno.h"
#include"sflow_diffusion_void.h"
#include"sflow_hydrostatic.h"
#include"sflow_eta.h"
#include"sflow_forcing.h"
#include"sflow_momentum_RK3.h"
#include"sflow_pjm_lin.h"
#include"sflow_ediff.h"
#include"sflow_amr_ship.h"
#include"6DOF_sflow.h"
#include"reefmg_core.h"
#include"reefmg2D.h"
#include"vec2D.h"
#include"ioflow_void.h"
#include"mgcslice1.h"
#include"mgcslice2.h"
#include"mgcslice4.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<cstdio>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
// floor(a/2^l) for negative a as well
inline int fsh(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }

inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b)?a:b;
}

// index of cell (i,j) in the lexer 2D arrays (wet, flagslice4)
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

inline bool inarr(const lexer *q, int ii, int jj)
{
    return ii>=q->imin && ii<q->imin+q->imax && jj>=q->jmin && jj<q->jmin+q->jmax;
}

// switches the MPI exchange of the ghostcell class off while patch kernels run
struct comms_off
{
    ghostcell *g;
    bool old;
    comms_off(ghostcell *gg) : g(gg) { old = g->set_comms(false); }
    ~comms_off() { g->set_comms(old); }
};
}

sflow_amr::sflow_amr(lexer *p, fdm2D *b, ghostcell *pgc, patchBC_interface *ppBC, sixdof *pp6dof,
                     sflow_HLL *pphll, sflow_momentum_RK3 *ppmom) : eps(1.0e-6)
{
    p0 = p;
    b0 = b;
    pgc0 = pgc;
    pBC = ppBC;
    p6dof = pp6dof;
    phll0 = pphll;
    pmom0 = ppmom;

    maxlev = p->A270;
    regrid_int = p->A271;
    nbuf = MAX(p->A272,0);
    tol_eta = p->A273;
    shore = p->A274;
    tile = MAX(p->A275,4);
    tile += tile%2;
    nest = 3;

    // cells computed beyond the patch box: the wet-dry step of a stage needs the new state
    // two cells further out, the discharge limiter (B 60) three
    // (non-hydrostatic without the shallow-water switch A 221: the deep flag needs the
    // wet state three cells further out)
    EXT = (p->B60>=1) ? 3 : 2;
    if(p->A220==1 && p->A221==0)
    EXT = 4;

    nh = (p->A220==1) ? 1 : 0;
    nh_it_total = nh_solves = 0;
    nh_it_last = 0;
    nh_rebuild0 = true;

    // moving body (X 10 2: direct forcing, X 10 3: pressure of the ship)
    shipmode = 0;
    ship6 = nullptr;
    if(p->X10==2 || p->X10==3)
    {
        ship6 = dynamic_cast<sixdof_sflow*>(pp6dof);
        if(ship6!=nullptr)
        shipmode = p->X10;
    }

    // boxes only: static patches
    if(tol_eta<=0.0 && shore==0 && (p->A278==0 || shipmode==0))
    regrid_int = 0;

    patches_total = 0;
    cells_total = 0;
    printcount_amr = 0;
    printtime_amr = 0.0;
    regrids = 0;
    m0 = 0.0;
    last_owner = 0;
    for(int k=0; k<9; ++k)
    tm[k]=0.0;

    match.resize(1);
    rmatch.resize(1);

    pflow_void = new ioflow_v(p,pgc,pBC);

    phll0->amr = this;
    phll0->amr_id = 0;
}

sflow_amr::~sflow_amr()
{
    for(auto c : P)
    free_patch(c);

    for(auto v : nh0_v)
    delete v;
    delete nhmg0;
    delete nhr0;
}

// --------------------------------------------------------------------- grids
sflow_amr::gh sflow_amr::grid(int g)
{
    gh r;
    if(g<0)
    {
        r.q = p0; r.b = b0; r.m = pmom0;
        r.oi = O0i; r.oj = O0j;
        return r;
    }
    sflow_amr_patch *c = P[g];
    r.q = c->pp; r.b = c->b; r.m = c->pmom;
    r.oi = c->I0-EXT; r.oj = c->J0-EXT;
    return r;
}

int sflow_amr::patch_at(int l, int I, int J)
{
    if(I<bxlo(l) || I>bxhi(l) || J<bylo(l) || J>byhi(l))
    return -3;

    if(l==0)
    return -1;

    int ti = I/tile - tti0[l];
    int tj = J/tile - ttj0[l];
    return tmap[l][ti*tny[l]+tj];
}

int sflow_amr::owner(int I, int J)
{
    if(I<0 || I>=GNX || J<0 || J>=GNY)
    return -1;

    int r = last_owner;
    if(I>=rbx0[r] && I<=rbx1[r] && J>=rby0[r] && J<=rby1[r])
    return r;

    for(r=0; r<(int)rbx0.size(); ++r)
    if(I>=rbx0[r] && I<=rbx1[r] && J>=rby0[r] && J<=rby1[r])
    {
        last_owner = r;
        return r;
    }

    return -1;
}

// level-0 solid flags of the rank box and a halo of FH cells (the ghostcell halo
// of flagslice4 is only one cell deep)
int sflow_amr::flag0(int I, int J)
{
    if(I<0 || I>=GNX || J<0 || J>=GNY)
    return -10;

    int ii = I-O0i+FH, jj = J-O0j+FH;
    if(ii<0 || jj<0 || ii>=NX0+2*FH || jj>=NY0+2*FH)
    return -10;

    return fl0[(size_t)ii*(NY0+2*FH)+jj];
}

void sflow_amr::build_flags(lexer *p)
{
    const int me = p->mpirank;
    const int ny = NY0+2*FH;
    fl0.assign((size_t)(NX0+2*FH)*ny,-10);

    vector<vector<int>> req(p->mpi_size), srv;
    vector<vector<int>> pos(p->mpi_size);

    for(int ii=0; ii<NX0+2*FH; ++ii)
    for(int jj=0; jj<ny; ++jj)
    {
        int I = ii-FH+O0i, J = jj-FH+O0j;
        int o = owner(I,J);
        if(o<0)
        continue;

        if(o==me)
        {
            fl0[(size_t)ii*ny+jj] = p->flagslice4[lij(p,I-O0i,J-O0j)];
            continue;
        }

        req[o].push_back(I); req[o].push_back(J);
        pos[o].push_back(ii*ny+jj);
    }

    xsetup(req,srv,2);

    sflow_amr_xplan X;
    for(int r=0; r<p->mpi_size; ++r)
    {
        if(!srv[r].empty())
        {
            X.speer.push_back(r);
            vector<double> v;
            for(size_t k=0; k<srv[r].size(); k+=2)
            v.push_back(double(p->flagslice4[lij(p,srv[r][k]-O0i,srv[r][k+1]-O0j)]));
            X.sbuf.push_back(v);
        }
        if(!req[r].empty())
        {
            X.rpeer.push_back(r);
            X.rcount.push_back((int)req[r].size()/2);
        }
    }
    xrun(X,1,7002);

    for(size_t k=0; k<X.rpeer.size(); ++k)
    {
        int r = X.rpeer[k];
        for(size_t m=0; m<pos[r].size(); ++m)
        fl0[pos[r][m]] = (int)X.rbuf[k][m];
    }
}

// --------------------------------------------------------------------- setup
void sflow_amr::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // scope
    int ok=1;
    if(p->A220!=0 && p->A220!=1) ok=0;
    if(p->A210==2) ok=0;
    if(p->A260!=0) ok=0;
    if(p->S10!=0) ok=0;
    if(p->X10>3) ok=0;
    if(p->W90!=0) ok=0;
    if(p->j_dir!=1) ok=0;

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"SFLOW AMR (A 270): only for A 220 0/1, A 210 3, A 260 0, S 10 0, X 10 0-3, W 90 0 and 2D grids -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    O0i = p->origin_i;
    O0j = p->origin_j;
    NX0 = p->knox;
    NY0 = p->knoy;
    GNX = p->gknox;
    GNY = p->gknoy;

    // rank boxes
    int mine[4] = {O0i, O0i+NX0-1, O0j, O0j+NY0-1};
    vector<int> all(4*p->mpi_size);
    MPI_Allgather(mine,4,MPI_INT,&all[0],4,MPI_INT,MPI_COMM_WORLD);
    rbx0.resize(p->mpi_size); rbx1.resize(p->mpi_size);
    rby0.resize(p->mpi_size); rby1.resize(p->mpi_size);
    for(int r=0; r<p->mpi_size; ++r)
    {
        rbx0[r]=all[4*r]; rbx1[r]=all[4*r+1];
        rby0[r]=all[4*r+2]; rby1[r]=all[4*r+3];
    }

    lev.assign(maxlev+1,vector<int>());
    gtile.assign(maxlev+1,vector<unsigned char>());
    gtnx.assign(maxlev+1,0);
    gtny.assign(maxlev+1,0);
    tmap.assign(maxlev+1,vector<int>());
    tti0.assign(maxlev+1,0); ttj0.assign(maxlev+1,0);
    tnx.assign(maxlev+1,0); tny.assign(maxlev+1,0);
    gplan.assign(maxlev+1,sflow_amr_xplan());
    gserve.assign(maxlev+1,vector<sflow_amr_fill>());
    grecv.assign(maxlev+1,vector<sflow_amr_fill*>());
    fplan.assign(maxlev+1,sflow_amr_xplan());
    fsend.assign(maxlev+1,vector<int>());
    frecv.assign(maxlev+1,vector<int>());

    for(int l=1; l<=maxlev; ++l)
    {
        gtnx[l] = ((GNX<<l)+tile-1)/tile;
        gtny[l] = ((GNY<<l)+tile-1)/tile;
    }
    build_tiles();
    build_flags(p);

    // bed and still water depth in the level-0 halo (the fine bed is prolonged from it)
    pgc->gcslparax(p,b->bed,4);
    pgc->gcslparax(p,b->depth,4);

    // no refinement next to in- and outflow boundaries and in the A 277 boxes
    forbid0.assign((size_t)GNX*GNY,0);
    auto forbid_around = [&](int ii, int jj)
    {
        int I=ii+O0i, J=jj+O0j;
        for(int a=MAX(I-4,0); a<=MIN(I+4,GNX-1); ++a)
        for(int d=MAX(J-4,0); d<=MIN(J+4,GNY-1); ++d)
        forbid0[(size_t)a*GNY+d]=1;
    };
    for(n=0; n<p->gcslin_count; ++n)
    forbid_around(p->gcslin[n][0],p->gcslin[n][1]);
    for(n=0; n<p->gcslout_count; ++n)
    forbid_around(p->gcslout[n][0],p->gcslout[n][1]);

    for(int q=0; q<p->A277; ++q)
    SLICELOOP4
    if(p->XP[IP]>=p->A277_xs[q] && p->XP[IP]<=p->A277_xe[q] && p->YP[JP]>=p->A277_ys[q] && p->YP[JP]<=p->A277_ye[q])
    forbid0[(size_t)(i+O0i)*GNY+(j+O0j)]=1;

    global_or(forbid0);

    mkdir("./REEF3D_SFLOW_AMR",0777);

    if(p->mpirank==0)
    {
        logout.open("./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_log.dat");
        logout<<"# count \t simtime \t dt \t patches \t cells \t mass \t rel. mass change"<<endl;
    }

    // initial hierarchy: every pass can add one level
    for(int it=0; it<maxlev; ++it)
    regrid(p,b,pgc,true);

    m0 = mass(p,b,pgc);

    if(p->mpirank==0)
    {
        cout<<"SFLOW AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined cells";
        if(regrid_int>0)
        cout<<", regrid every "<<regrid_int<<" steps";
        cout<<endl;
        if(shipmode>0)
        cout<<"SFLOW AMR: moving body X 10 "<<shipmode<<" on the patches"<<(p->A278>0?", refinement around the body (A 278)":"")<<endl;
        if(p->A278>0 && shipmode==0)
        cout<<"SFLOW AMR: A 278 needs a moving body (X 10 2/3) -- ignored"<<endl;
    }
}

// allocate the patch and its SFLOW objects; the bed and the state are set later
sflow_amr_patch* sflow_amr::make_patch(lexer *p, ghostcell *pgc, int l, int I0, int I1, int J0, int J1)
{
    sflow_amr_patch *c = new sflow_amr_patch;
    c->lev = l;
    c->I0=I0; c->I1=I1; c->J0=J0; c->J1=J1;
    c->nx = I1-I0+1;
    c->ny = J1-J0+1;
    c->fresh = true;

    c->pp = new lexer(*p,1);
    build_lexer(p,*c);

    {
    comms_off guard(pgc);

    build_bc(pgc,*c);

    lexer *pp = c->pp;
    c->b = new fdm2D(pp);
    c->phll = new sflow_HLL(pp,pgc,pBC);
    c->phll->amr = this;
    c->pss = new sflow_signal_speed(pp);

    if(p->A211<=3)
    c->precon = new sflow_reconstruct_hires(pp,pBC);
    if(p->A211>=4)
    c->precon = new sflow_reconstruct_weno(pp,pBC,1);

    // diffusion: explicit on the patches (A 212 2 is implicit on level 0 only)
    if(p->A212>=1)
    c->pdiff = new sflow_ediff(pp);
    else
    c->pdiff = new sflow_diffusion_void(pp);
    if(nh==1)
    {
    c->pnh = new sflow_pjm_lin(pp,c->b,pBC);
    c->ppress = c->pnh;
    }
    else
    c->ppress = new sflow_hydrostatic(pp,c->b,pBC);
    c->pfsf = new sflow_eta(pp,c->b,pgc,pBC);
    c->psfdf = new sflow_forcing(pp);
    sixdof *p6 = p6dof;
    if(shipmode>0)
    {
    c->pship = new sflow_amr_ship(pp);
    p6 = c->pship;
    }
    c->pmom = new sflow_momentum_RK3(pp,c->b,pgc,c->phll,c->pss,c->precon,c->pdiff,c->ppress,nullptr,nullptr,pflow_void,c->pfsf,c->psfdf,p6);
    c->pmom->nh_defer = (nh==1);
    }

    for(int ip=0; ip<5; ++ip)
    {
        c->rec[ip][0].assign(c->ny,0.0);
        c->rec[ip][1].assign(c->ny,0.0);
        c->rec[ip][2].assign(c->nx,0.0);
        c->rec[ip][3].assign(c->nx,0.0);
    }

    return c;
}

void sflow_amr::free_patch(sflow_amr_patch *c)
{
    lexer *pp = c->pp;

    for(auto v : c->nv)
    delete v;
    delete c->mg;

    delete c->pmom;
    delete c->pship;
    delete c->psfdf;
    delete c->pfsf;
    delete c->ppress;
    delete c->pdiff;
    delete c->precon;
    delete c->pss;
    delete c->phll;
    delete c->b;

    const int nxa = pp->knox+1+4*marge;
    const int nya = pp->knoy+1+4*marge;
    const int nsl = pp->imax*pp->jmax;

    pp->del_Darray(pp->XN,nxa); pp->del_Darray(pp->XP,nxa); pp->del_Darray(pp->DXN,nxa); pp->del_Darray(pp->DXP,nxa);
    pp->del_Darray(pp->YN,nya); pp->del_Darray(pp->YP,nya); pp->del_Darray(pp->DYN,nya); pp->del_Darray(pp->DYP,nya);
    pp->del_Iarray(pp->flagslice1,nsl);
    pp->del_Iarray(pp->flagslice2,nsl);
    pp->del_Iarray(pp->flagslice4,nsl);
    pp->del_Iarray(pp->wet,nsl);
    pp->del_Iarray(pp->wet_n,nsl);
    pp->del_Iarray(pp->deep,nsl);
    pp->del_Iarray(pp->IOSL,nsl);
    pp->del_Iarray(pp->sizeS1,5);
    pp->del_Iarray(pp->sizeS2,5);
    pp->del_Iarray(pp->sizeS4,5);
    pp->del_Iarray(pp->gcbsl1,pp->gcbsl1_count,6);
    pp->del_Iarray(pp->gcbsl2,pp->gcbsl2_count,6);
    pp->del_Iarray(pp->gcbsl4,pp->gcbsl4_count,6);
    pp->del_Iarray(pp->dgcsl1,pp->dgcsl1_count,3);
    pp->del_Iarray(pp->dgcsl2,pp->dgcsl2_count,3);
    pp->del_Iarray(pp->dgcsl4,pp->dgcsl4_count,3);

    delete pp;
    delete c;
}

// patch geometry: nodes of the global level-l grid, subdividing the level-0 nodes
void sflow_amr::build_lexer(lexer *p, sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    const int ms = pp->margin;   // ghost layers of the slices
    const int m = marge;         // offset of the coordinate arrays (IP = i+marge)
    const int l = c.lev;
    const double rl = double(1<<l);

    pp->knox = c.nx+2*EXT;
    pp->knoy = c.ny+2*EXT;
    pp->knoz = p->knoz;
    pp->imin = -ms;
    pp->jmin = -ms;
    pp->imax = pp->knox+2*ms;
    pp->jmax = pp->knoy+2*ms;
    pp->kmin = p->kmin;
    pp->kmax = p->kmax;
    pp->kmaxF = p->kmaxF;

    pp->i_dir = p->i_dir; pp->j_dir = p->j_dir; pp->k_dir = p->k_dir;
    pp->x_dir = p->x_dir; pp->y_dir = p->y_dir; pp->z_dir = p->z_dir;

    pp->ulast = pp->vlast = pp->wlast = pp->flast = 0;
    pp->ulastsflow = 0;

    // coordinates
    const int nxa = pp->knox+1+4*m;
    const int nya = pp->knoy+1+4*m;
    const int pnxa = p->knox+1+4*m;
    const int pnya = p->knoy+1+4*m;

    pp->Darray(pp->XN,nxa); pp->Darray(pp->XP,nxa); pp->Darray(pp->DXN,nxa); pp->Darray(pp->DXP,nxa);
    pp->Darray(pp->YN,nya); pp->Darray(pp->YP,nya); pp->Darray(pp->DYN,nya); pp->Darray(pp->DYP,nya);
    pp->ZN = p->ZN; pp->ZP = p->ZP; pp->DZN = p->DZN; pp->DZP = p->DZP;

    auto x0 = [&](int K)   // global level-0 node K, linear extrapolation outside the rank arrays
    {
        int a = K-O0i+m;
        if(a<0)      return p->XN[0] + a*(p->XN[1]-p->XN[0]);
        if(a>=pnxa)  return p->XN[pnxa-1] + (a-pnxa+1)*(p->XN[pnxa-1]-p->XN[pnxa-2]);
        return p->XN[a];
    };
    auto y0 = [&](int K)
    {
        int a = K-O0j+m;
        if(a<0)      return p->YN[0] + a*(p->YN[1]-p->YN[0]);
        if(a>=pnya)  return p->YN[pnya-1] + (a-pnya+1)*(p->YN[pnya-1]-p->YN[pnya-2]);
        return p->YN[a];
    };

    for(int a=0; a<nxa; ++a)
    {
        int N = a-m-EXT+c.I0;      // global level-l node
        int K = fsh(N,l);
        int r = N-(K<<l);
        pp->XN[a] = (r==0) ? x0(K) : x0(K) + (x0(K+1)-x0(K))*(double(r)/rl);
    }
    for(int a=0; a<nya; ++a)
    {
        int N = a-m-EXT+c.J0;
        int K = fsh(N,l);
        int r = N-(K<<l);
        pp->YN[a] = (r==0) ? y0(K) : y0(K) + (y0(K+1)-y0(K))*(double(r)/rl);
    }

    for(int a=0; a<nxa-1; ++a) { pp->XP[a] = 0.5*(pp->XN[a]+pp->XN[a+1]); pp->DXN[a] = pp->XN[a+1]-pp->XN[a]; }
    for(int a=0; a<nya-1; ++a) { pp->YP[a] = 0.5*(pp->YN[a]+pp->YN[a+1]); pp->DYN[a] = pp->YN[a+1]-pp->YN[a]; }
    pp->XP[nxa-1] = pp->XP[nxa-2] + pp->DXN[nxa-2]; pp->DXN[nxa-1] = pp->DXN[nxa-2];
    pp->YP[nya-1] = pp->YP[nya-2] + pp->DYN[nya-2]; pp->DYN[nya-1] = pp->DYN[nya-2];
    for(int a=0; a<nxa-1; ++a) pp->DXP[a] = pp->XP[a+1]-pp->XP[a];
    for(int a=0; a<nya-1; ++a) pp->DYP[a] = pp->YP[a+1]-pp->YP[a];
    pp->DXP[nxa-1] = pp->DXP[nxa-2];
    pp->DYP[nya-1] = pp->DYP[nya-2];

    // mean spacing: exact fraction of level 0, independent of the patch size
    pp->DXM = p->DXM/rl;
    pp->DYM = pp->DXM;
    pp->DXD = p->DXD;
    pp->DYD = p->DYD;

    // flags from level 0, wet arrays
    const int nsl = pp->imax*pp->jmax;
    pp->Iarray(pp->flagslice1,nsl);
    pp->Iarray(pp->flagslice2,nsl);
    pp->Iarray(pp->flagslice4,nsl);
    pp->Iarray(pp->wet,nsl);
    pp->Iarray(pp->wet_n,nsl);
    pp->Iarray(pp->deep,nsl);

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        int I = ii-EXT+c.I0, J = jj-EXT+c.J0;
        pp->flagslice4[lij(pp,ii,jj)] = flag0(fsh(I,l),fsh(J,l));
    }

    // face flags as mgcslice1/2::makemgc, over the whole array
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        int q = lij(pp,ii,jj);
        int f4 = pp->flagslice4[q];
        pp->flagslice1[q] = f4;
        pp->flagslice2[q] = f4;

        if(f4>0 && inarr(pp,ii+1,jj) && pp->flagslice4[lij(pp,ii+1,jj)]<0)
        pp->flagslice1[q] = pp->flagslice4[lij(pp,ii+1,jj)];

        if(f4>0 && inarr(pp,ii,jj+1) && pp->flagslice4[lij(pp,ii,jj+1)]<0)
        pp->flagslice2[q] = pp->flagslice4[lij(pp,ii,jj+1)];
    }

    for(int q=0; q<nsl; ++q)
    {
        pp->wet[q] = pp->wet_n[q] = pp->deep[q] = 0;
    }

    // boundary lists are built in build_bc; no partition neighbours
    pp->gcbsl1 = pp->gcbsl2 = pp->gcbsl3 = pp->gcbsl4 = pp->gcbsl4a = nullptr;
    pp->dgcsl1 = pp->dgcsl2 = pp->dgcsl3 = pp->dgcsl4 = nullptr;
    pp->gcslin = pp->gcslout = pp->gcslawa1 = pp->gcslawa2 = nullptr;
    pp->gcbsl1_count = pp->gcbsl2_count = pp->gcbsl3_count = pp->gcbsl4_count = pp->gcbsl4a_count = 0;
    pp->gcslin_count = pp->gcslout_count = 0;
    pp->gcslawa1_count = pp->gcslawa2_count = 0;
    pp->dgcsl1_count = pp->dgcsl2_count = pp->dgcsl3_count = pp->dgcsl4_count = 0;
    pp->gcsldfeta4_count = pp->gcsldfbed4_count = 0;
    pp->gcslpara1_count = pp->gcslpara2_count = pp->gcslpara3_count = pp->gcslpara4_count = 0;
    pp->gcslparaco1_count = pp->gcslparaco2_count = pp->gcslparaco3_count = pp->gcslparaco4_count = 0;
    pp->nb1 = pp->nb2 = pp->nb3 = pp->nb4 = pp->nb5 = pp->nb6 = -2;
    pp->mpi_size = 1;
    pp->mpirank = p->mpirank;

    // the domain boundary never lies on a patch: all boundary faces are walls (bc 21)
    pp->origin_i = pp->origin_j = pp->origin_k = -100000000;
    pp->gknox = p->gknox; pp->gknoy = p->gknoy; pp->gknoz = p->gknoz;

    pp->vec2Dlength = nsl;
    pp->veclength = 0;
    pp->cellnum2D = pp->knox*pp->knoy;
    pp->slicenum = nsl;

    pp->global_xmin = p->global_xmin; pp->global_ymin = p->global_ymin;
    pp->global_xmax = p->global_xmax; pp->global_ymax = p->global_ymax;
    pp->originx = p->originx; pp->originy = p->originy;

    pp->wd = p->wd;
    pp->phimean = p->phimean;
    pp->phiout = p->phiout;
    pp->dt = p->dt;
    pp->dt_old = p->dt_old;
    pp->simtime = p->simtime;
    pp->count = p->count;
}

// wall lists of the patch, as driver::makegrid2D for level 0
void sflow_amr::build_bc(ghostcell *pgc, sflow_amr_patch &c)
{
    lexer *pp = c.pp;

    pp->gridini2D();

    mgcslice1 m1(pp);
    mgcslice2 m2(pp);
    mgcslice4 m4(pp);

    m1.gcb_seed(pp);
    m2.gcb_seed(pp);
    m4.gcb_seed(pp);

    pgc->gcsl_setbc1(pp);
    pgc->gcsl_setbc2(pp);
    pgc->gcsl_setbc4(pp);

    pgc->gcsl_setbcio(pp);

    pgc->dgcslini1(pp);
    pgc->dgcslini2(pp);
    pgc->dgcslini4(pp);
}

// fine bed: limited-linear prolongation, level by level from level 0 (children average = parent);
// solid cells and cells next to solids take the parent value without slopes
double sflow_amr::bed_at(int l, int I, int J)
{
    if(l==0)
    {
        // cells outside the domain take the bed of the nearest domain cell
        I = MAX(MIN(I,GNX-1),0);
        J = MAX(MIN(J,GNY-1),0);
        int ii = MAX(MIN(I-O0i, p0->imin+p0->imax-1), p0->imin);
        int jj = MAX(MIN(J-O0j, p0->jmin+p0->jmax-1), p0->jmin);
        return b0->bed(ii,jj);
    }

    const uint64_t key = (uint64_t(l)<<58) ^ (uint64_t(uint32_t(I+(1<<28)))<<29) ^ uint64_t(uint32_t(J+(1<<28)));
    auto it = bedmemo.find(key);
    if(it!=bedmemo.end())
    return it->second;

    int Ic = fsh(I,1), Jc = fsh(J,1);
    double bc = bed_at(l-1,Ic,Jc);

    if(flag0(fsh(I,l),fsh(J,l))<0)
    {
        bedmemo[key]=bc;
        return bc;
    }

    const int lc = l-1;
    auto fl = [&](int a, int d) { return flag0(fsh(a,lc),fsh(d,lc)); };

    double sx=0.0, sy=0.0;
    if(fl(Ic+1,Jc)>0 && fl(Ic-1,Jc)>0)
    sx = mmod(bed_at(lc,Ic+1,Jc)-bc, bc-bed_at(lc,Ic-1,Jc));
    if(fl(Ic,Jc+1)>0 && fl(Ic,Jc-1)>0)
    sy = mmod(bed_at(lc,Ic,Jc+1)-bc, bc-bed_at(lc,Ic,Jc-1));

    double ox = (I-2*Ic==0) ? -0.25 : 0.25;
    double oy = (J-2*Jc==0) ? -0.25 : 0.25;
    double r = bc + sx*ox + sy*oy;
    bedmemo[key]=r;
    return r;
}

void sflow_amr::build_tiles()
{
    for(int l=1; l<=maxlev; ++l)
    {
        tti0[l] = bxlo(l)/tile;
        ttj0[l] = bylo(l)/tile;
        tnx[l] = bxhi(l)/tile - tti0[l] + 1;
        tny[l] = byhi(l)/tile - ttj0[l] + 1;
        tmap[l].assign(tnx[l]*tny[l],-2);

        for(int n : lev[l])
        {
            sflow_amr_patch *c = P[n];
            for(int ti=c->I0/tile; ti<=c->I1/tile; ++ti)
            for(int tj=c->J0/tile; tj<=c->J1/tile; ++tj)
            tmap[l][(ti-tti0[l])*tny[l]+(tj-ttj0[l])] = n;
        }
    }
}

void sflow_amr::global_or(vector<unsigned char> &v)
{
    if(p0->mpi_size>1 && !v.empty())
    MPI_Allreduce(MPI_IN_PLACE,&v[0],(int)v.size(),MPI_UNSIGNED_CHAR,MPI_BOR,MPI_COMM_WORLD);
}

// --------------------------------------------------------------------- refinement flags
// cells of level l on this rank that need level l+1 (evaluated on the level-l grids)
void sflow_amr::tag_level(int l, vector<unsigned char> &M)
{
    const int lf = l+1;
    const int T = tile;
    const double hfilm = 3.0*p0->A244;

    auto mark = [&](int I, int J, int buf)
    {
        int a0 = MAX((2*(I-buf))/T,0), a1 = MIN((2*(I+buf)+1)/T,gtnx[lf]-1);
        int d0 = MAX((2*(J-buf))/T,0), d1 = MIN((2*(J+buf)+1)/T,gtny[lf]-1);
        if(2*(I-buf)<0) a0=0;
        if(2*(J-buf)<0) d0=0;
        for(int a=a0; a<=a1; ++a)
        for(int d=d0; d<=d1; ++d)
        M[(size_t)a*gtny[lf]+d]=1;
    };

    auto test = [&](lexer *q, fdm2D *bb, slice &WL, int ii, int jj)
    {
        if(q->flagslice4[lij(q,ii,jj)]<0)
        return false;

        int wc = q->wet[lij(q,ii,jj)];
        double ec = WL(ii,jj)-bb->depth(ii,jj);
        const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};

        for(int k=0; k<4; ++k)
        {
            int a=ii+di[k], d=jj+dj[k];
            if(q->flagslice4[lij(q,a,d)]<0)
            continue;

            int wn = q->wet[lij(q,a,d)];

            if(tol_eta>0.0 && wc==1 && wn==1)
            if(fabs(ec-(WL(a,d)-bb->depth(a,d)))>tol_eta)
            return true;

            if(shore>0 && wc==1 && wn==0 && WL(ii,jj)>hfilm)
            return true;
        }
        return false;
    };

    if(l==0)
    {
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(test(p0,b0,b0->WL,ii,jj))
        mark(ii+O0i,jj+O0j,nbuf);
        return;
    }

    for(int n : lev[l])
    {
        sflow_amr_patch *c = P[n];
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(test(c->pp,c->b,c->b->WL,ii,jj))
        mark(ii-EXT+c->I0,jj-EXT+c->J0,nbuf);
    }
}

// --------------------------------------------------------------------- regridding
void sflow_amr::regrid(lexer *p, fdm2D *b, ghostcell *pgc, bool initial)
{
    const int T = tile;
    const int me = p->mpirank;
    const int old_total = patches_total;

    // current hierarchy: final state with filled ghost cells
    cache_stage(0);
    for(int l=1; l<=maxlev; ++l)
    fill_level(pgc,l,0);

    // ---- tile maps, finest first: flags of the next coarser level, the A 276 boxes
    //      and the footprint of the finer level grown by the nesting width
    vector<vector<unsigned char>> M(maxlev+1);

    if(shipmode>0 && p->A278>0)
    ship_setup(p);

    auto tile_forbidden = [&](int l, int ti, int tj)
    {
        int a0 = (ti*T)>>l, a1 = (MIN((ti+1)*T,GNX<<l)-1)>>l;
        int d0 = (tj*T)>>l, d1 = (MIN((tj+1)*T,GNY<<l)-1)>>l;
        for(int a=a0; a<=a1; ++a)
        for(int d=d0; d<=d1; ++d)
        if(forbid0[(size_t)a*GNY+d])
        return true;
        return false;
    };

    for(int l=maxlev; l>=1; --l)
    {
        M[l].assign((size_t)gtnx[l]*gtny[l],0);

        tag_level(l-1,M[l]);

        for(int q=0; q<p->A276; ++q)
        SLICELOOP4
        if(p->XP[IP]>=p->A276_xs[q] && p->XP[IP]<=p->A276_xe[q] && p->YP[JP]>=p->A276_ys[q] && p->YP[JP]<=p->A276_ye[q])
        {
            int I=i+O0i, J=j+O0j;
            for(int ti=(I<<l)/T; ti<=(((I+1)<<l)-1)/T; ++ti)
            for(int tj=(J<<l)/T; tj<=(((J+1)<<l)-1)/T; ++tj)
            M[l][(size_t)ti*gtny[l]+tj]=1;
        }

        if(shipmode>0 && p->A278>0)
        SLICELOOP4
        if(ship_zone(p->XP[IP],p->YP[JP]))
        {
            int I=i+O0i, J=j+O0j;
            for(int ti=(I<<l)/T; ti<=(((I+1)<<l)-1)/T; ++ti)
            for(int tj=(J<<l)/T; tj<=(((J+1)<<l)-1)/T; ++tj)
            M[l][(size_t)ti*gtny[l]+tj]=1;
        }

        if(l<maxlev)
        for(int ti=0; ti<gtnx[l+1]; ++ti)
        for(int tj=0; tj<gtny[l+1]; ++tj)
        if(M[l+1][(size_t)ti*gtny[l+1]+tj])
        {
            int a0 = (ti*T)/2-nest, a1 = ((ti+1)*T-1)/2+nest;
            int d0 = (tj*T)/2-nest, d1 = ((tj+1)*T-1)/2+nest;
            a0 = MAX(a0,0)/T; a1 = MIN(a1,(GNX<<l)-1)/T;
            d0 = MAX(d0,0)/T; d1 = MIN(d1,(GNY<<l)-1)/T;
            for(int a=a0; a<=a1; ++a)
            for(int d=d0; d<=d1; ++d)
            M[l][(size_t)a*gtny[l]+d]=1;
        }

        global_or(M[l]);

        for(int ti=0; ti<gtnx[l]; ++ti)
        for(int tj=0; tj<gtny[l]; ++tj)
        if(M[l][(size_t)ti*gtny[l]+tj] && tile_forbidden(l,ti,tj))
        M[l][(size_t)ti*gtny[l]+tj]=0;
    }

    // proper nesting (the forbidden tiles can break it): coarse first
    for(int l=2; l<=maxlev; ++l)
    for(int ti=0; ti<gtnx[l]; ++ti)
    for(int tj=0; tj<gtny[l]; ++tj)
    if(M[l][(size_t)ti*gtny[l]+tj])
    {
        int a0 = (ti*T)/2-nest, a1 = ((ti+1)*T-1)/2+nest;
        int d0 = (tj*T)/2-nest, d1 = ((tj+1)*T-1)/2+nest;
        a0 = MAX(a0,0)/T; a1 = MIN(a1,(GNX<<(l-1))-1)/T;
        d0 = MAX(d0,0)/T; d1 = MIN(d1,(GNY<<(l-1))-1)/T;
        bool ok=true;
        for(int a=a0; a<=a1 && ok; ++a)
        for(int d=d0; d<=d1 && ok; ++d)
        if(!M[l-1][(size_t)a*gtny[l-1]+d])
        ok=false;

        if(!ok)
        M[l][(size_t)ti*gtny[l]+tj]=0;
    }

    // ---- patches: marked tiles merged into rectangles, cut at the rank box
    vector<sflow_amr_patch*> oldP = P;
    vector<char> keep(oldP.size(),0);
    vector<sflow_amr_patch*> newP;
    vector<vector<int>> newlev(maxlev+1);

    for(int l=1; l<=maxlev; ++l)
    {
        gtile[l] = M[l];

        struct R { int ta,tb,tj0,tj1; };
        vector<R> open, done;
        int t0i=bxlo(l)/T, t1i=bxhi(l)/T, t0j=bylo(l)/T, t1j=byhi(l)/T;

        for(int tj=t0j; tj<=t1j; ++tj)
        {
            vector<R> next;
            int ti=t0i;
            while(ti<=t1i)
            {
                if(!M[l][(size_t)ti*gtny[l]+tj]) { ++ti; continue; }
                int ta=ti;
                while(ti<=t1i && M[l][(size_t)ti*gtny[l]+tj]) ++ti;
                int tb=ti-1;

                bool ext=false;
                for(auto &o : open)
                if(o.ta==ta && o.tb==tb && o.tj1==tj-1)
                {
                    o.tj1=tj;
                    next.push_back(o);
                    o.ta=-1;
                    ext=true;
                    break;
                }
                if(!ext)
                next.push_back({ta,tb,tj,tj});
            }
            for(auto &o : open)
            if(o.ta>=0)
            done.push_back(o);
            open = next;
        }
        for(auto &o : open)
        done.push_back(o);

        for(auto &o : done)
        {
            int I0 = MAX(o.ta*T,bxlo(l)), I1 = MIN((o.tb+1)*T-1,bxhi(l));
            int J0 = MAX(o.tj0*T,bylo(l)), J1 = MIN((o.tj1+1)*T-1,byhi(l));

            // only solid cells: no patch
            int fluid=0;
            for(int a=(I0>>l); a<=(I1>>l) && fluid==0; ++a)
            for(int d=(J0>>l); d<=(J1>>l) && fluid==0; ++d)
            if(flag0(a,d)>0)
            fluid=1;
            if(fluid==0)
            continue;

            sflow_amr_patch *c=nullptr;
            for(size_t k=0; k<oldP.size(); ++k)
            if(!keep[k] && oldP[k]->lev==l && oldP[k]->I0==I0 && oldP[k]->I1==I1 && oldP[k]->J0==J0 && oldP[k]->J1==J1)
            {
                keep[k]=1;
                c = oldP[k];
                c->fresh = false;
                break;
            }

            if(c==nullptr)
            c = make_patch(p,pgc,l,I0,I1,J0,J1);

            newlev[l].push_back((int)newP.size());
            newP.push_back(c);
        }
    }

    vector<sflow_amr_patch*> gone;
    for(size_t k=0; k<oldP.size(); ++k)
    if(!keep[k])
    gone.push_back(oldP[k]);

    P = newP;
    lev = newlev;
    for(int n=0; n<(int)P.size(); ++n)
    P[n]->phll->amr_id = n+1;

    build_tiles();

    // ---- bed of the new patches: prolonged from level 0, cells of other ranks from their owner
    {
    vector<vector<int>> req(p->mpi_size), srv;
    vector<vector<double*>> dst(p->mpi_size);

    for(auto c : P)
    if(c->fresh)
    {
        lexer *pp = c->pp;
        const int l = c->lev;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            int I = ii-EXT+c->I0, J = jj-EXT+c->J0;
            int o = owner(fsh(I,l),fsh(J,l));

            if(o<0 || o==me)
            c->b->bed(ii,jj) = bed_at(l,I,J);

            if(o>=0 && o!=me)
            {
                req[o].push_back(l); req[o].push_back(I); req[o].push_back(J);
                dst[o].push_back(&c->b->bed(ii,jj));
            }
        }
    }

    xsetup(req,srv,3);

    sflow_amr_xplan X;
    for(int r=0; r<p->mpi_size; ++r)
    {
        if(!srv[r].empty())
        {
            X.speer.push_back(r);
            vector<double> v;
            for(size_t k=0; k<srv[r].size(); k+=3)
            v.push_back(bed_at(srv[r][k],srv[r][k+1],srv[r][k+2]));
            X.sbuf.push_back(v);
        }
        if(!req[r].empty())
        {
            X.rpeer.push_back(r);
            X.rcount.push_back((int)req[r].size()/3);
        }
    }
    xrun(X,1,7001);

    for(size_t k=0; k<X.rpeer.size(); ++k)
    {
        int r = X.rpeer[k];
        for(size_t m=0; m<dst[r].size(); ++m)
        *dst[r][m] = X.rbuf[k][m];
    }

    for(auto c : P)
    if(c->fresh)
    {
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            pb->bed0(ii,jj) = pb->bed(ii,jj);
            pb->topobed(ii,jj) = pb->bed(ii,jj);
            pb->solidbed(ii,jj) = pb->bed(ii,jj);
            pb->depth(ii,jj) = p->wd - pb->bed(ii,jj);
            pb->ks(ii,jj) = p->B50;
        }

        // face depth as sflow_reconstruct::reconstruct_WL
        for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        pb->dfx(ii,jj) = 0.5*(pb->depth(ii+1,jj)+pb->depth(ii,jj));

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
        pb->dfy(ii,jj) = 0.5*(pb->depth(ii,jj+1)+pb->depth(ii,jj));
    }
    }

    // ---- fill, restriction and flux matching plans of the new hierarchy
    build_plans(pgc);
    cache_stage(0);

    // ---- state of the new patches, coarse to fine
    for(int l=1; l<=maxlev; ++l)
    {
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        if(P[n]->fresh)
        ini_patch_state(p,pgc,*P[n],oldP);
        }

        fill_level(pgc,l,0);
    }

    for(auto c : gone)
    free_patch(c);

    for(auto c : P)
    c->fresh = false;

    // level 0 consistent with the patches
    {
    comms_off guard(pgc);
    for(int l=maxlev; l>=1; --l)
    for(int n : lev[l])
    restrict_patch(p,*P[n],2);
    }

    long cells=0;
    for(auto c : P)
    cells += (long)c->nx*c->ny;

    patches_total = pgc->globalisum((int)P.size());
    nlevg.assign(maxlev+1,0);
    for(int l=1; l<=maxlev; ++l)
    nlevg[l] = pgc->globalisum((int)lev[l].size());
    cells_total = (long)pgc->globalsum(double(cells));

    if(patches_total>0 || old_total>0)
    exchange_level0(p,b,pgc,2);

    ++regrids;
}

// interior of a new patch: old patches of the same level where they existed, otherwise
// prolonged from the new coarser level (exact mass and momentum of fully wet parents)
void sflow_amr::ini_patch_state(lexer *p, ghostcell *pgc, sflow_amr_patch &c, vector<sflow_amr_patch*> &oldP)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;
    const int l = c.lev;
    const double wd = p->A244;
    double v[NV];

    // defaults for all cells (cells outside the interior are filled afterwards)
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        b->WL(ii,jj) = wd;
        b->UH(ii,jj) = b->VH(ii,jj) = b->WH(ii,jj) = 0.0;
        b->eta(ii,jj) = wd - b->depth(ii,jj);
        b->U(ii,jj) = b->V(ii,jj) = b->W(ii,jj) = 0.0;
        b->press(ii,jj) = 0.0;
        pp->wet[lij(pp,ii,jj)] = 0;
        pp->deep[lij(pp,ii,jj)] = 0;
    }

    for(int bi=0; bi<c.nx/2; ++bi)
    for(int bj=0; bj<c.ny/2; ++bj)
    {
        int i0 = EXT+2*bi, j0 = EXT+2*bj;
        int I = c.I0+2*bi, J = c.J0+2*bj;

        if(pp->flagslice4[lij(pp,i0,j0)]<0)
        continue;

        // old patch of the same level
        sflow_amr_patch *o=nullptr;
        for(auto q : oldP)
        if(q->lev==l && I>=q->I0 && I<=q->I1 && J>=q->J0 && J<=q->J1)
        {
            o=q;
            break;
        }

        if(o!=nullptr && o!=&c)
        {
            lexer *qp = o->pp;
            for(int a=0;a<2;++a)
            for(int d=0;d<2;++d)
            {
                int si = I+a-o->I0+EXT, sj = J+d-o->J0+EXT;
                int ii = i0+a, jj = j0+d;
                b->WL(ii,jj) = o->b->WL(si,sj);
                b->UH(ii,jj) = o->b->UH(si,sj);
                b->VH(ii,jj) = o->b->VH(si,sj);
                b->eta(ii,jj) = o->b->eta(si,sj);
                b->U(ii,jj) = o->b->U(si,sj);
                b->V(ii,jj) = o->b->V(si,sj);
                b->WH(ii,jj) = o->b->WH(si,sj);
                b->W(ii,jj) = o->b->W(si,sj);
                b->press(ii,jj) = o->b->press(si,sj);
                pp->wet[lij(pp,ii,jj)] = qp->wet[lij(qp,si,sj)];
                pp->deep[lij(pp,ii,jj)] = qp->deep[lij(qp,si,sj)];
            }
            continue;
        }

        // prolongation
        int Ic = I>>1, Jc = J>>1;
        int g = patch_at(l-1,Ic,Jc);
        if(g<-1)
        continue;

        gh G = grid(g);
        int ic = Ic-G.oi, jc = Jc-G.oj;

        for(int a=0;a<2;++a)
        for(int d=0;d<2;++d)
        {
            int ii = i0+a, jj = j0+d;
            prolong(g,ic,jc,a==0?-1:1,d==0?-1:1,b->depth(ii,jj),v);
            b->WL(ii,jj)=v[0]; b->UH(ii,jj)=v[1]; b->VH(ii,jj)=v[2]; b->WH(ii,jj)=v[8];
            b->press(ii,jj) = G.b->press(ic,jc);
            pp->wet[lij(pp,ii,jj)]=(int)v[6];
        }

        int nw = pp->wet[lij(pp,i0,j0)] + pp->wet[lij(pp,i0+1,j0)] + pp->wet[lij(pp,i0,j0+1)] + pp->wet[lij(pp,i0+1,j0+1)];

        slice &WLp = *sWL[g+1], &UHp = *sUH[g+1], &VHp = *sVH[g+1];

        if(nw==4 && G.q->wet[lij(G.q,ic,jc)]==1)
        {
            double sw = b->WL(i0,j0)+b->WL(i0+1,j0)+b->WL(i0,j0+1)+b->WL(i0+1,j0+1);
            double su = b->UH(i0,j0)+b->UH(i0+1,j0)+b->UH(i0,j0+1)+b->UH(i0+1,j0+1);
            double sv = b->VH(i0,j0)+b->VH(i0+1,j0)+b->VH(i0,j0+1)+b->VH(i0+1,j0+1);
            double sh = b->WH(i0,j0)+b->WH(i0+1,j0)+b->WH(i0,j0+1)+b->WH(i0+1,j0+1);
            double du = 4.0*UHp(ic,jc) - su;
            double dv = 4.0*VHp(ic,jc) - sv;
            double dw = 4.0*(*sWH[g+1])(ic,jc) - sh;
            double fw = 4.0*WLp(ic,jc)/sw;

            for(int a=0;a<2;++a)
            for(int d=0;d<2;++d)
            {
                b->UH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*du;
                b->VH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dv;
                b->WH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dw;
                b->WL(i0+a,j0+d) *= fw;
            }
        }

        for(int a=0;a<2;++a)
        for(int d=0;d<2;++d)
        {
            int ii = i0+a, jj = j0+d;
            int w = pp->wet[lij(pp,ii,jj)];
            double wlvl = fabs(b->WL(ii,jj))>wd ? b->WL(ii,jj) : 1.0e20;
            b->eta(ii,jj) = b->WL(ii,jj) - b->depth(ii,jj);
            b->U(ii,jj) = w==1 ? b->UH(ii,jj)/wlvl : 0.0;
            b->V(ii,jj) = w==1 ? b->VH(ii,jj)/wlvl*p->y_dir : 0.0;
            b->W(ii,jj) = w==1 ? b->WH(ii,jj)/wlvl : 0.0;
            pp->deep[lij(pp,ii,jj)] = w;
        }
    }

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        b->eta_n(ii,jj) = b->eta(ii,jj);
        b->hp(ii,jj) = b->WL(ii,jj);
        b->breaking(ii,jj) = 0;
        pp->wet_n[lij(pp,ii,jj)] = pp->wet[lij(pp,ii,jj)];
    }

    pp->dt = p->dt;
    pp->dt_old = p->dt_old;
    pp->simtime = p->simtime;
    pp->count = p->count;
}

// --------------------------------------------------------------------- plans
void sflow_amr::build_plans(ghostcell *pgc)
{
    const int me = p0->mpirank;
    const int np = p0->mpi_size;

    match.assign(P.size()+1,vector<sflow_amr_match>());
    rmatch.assign(P.size()+1,vector<sflow_amr_match>());

    // local source of a cell (l,I,J) owned by this rank: a patch of level l or the coarser level
    int violations=0;
    auto source = [&](int l, int I, int J, sflow_amr_fill &f)
    {
        int q = patch_at(l,I,J);
        if(q>=0)
        {
            f.kind = 0;
            f.g = q;
            f.si = I-P[q]->I0+EXT;
            f.sj = J-P[q]->J0+EXT;
            return true;
        }

        int Ic = I>>1, Jc = J>>1;
        int g = patch_at(l-1,Ic,Jc);
        if(g<-1)
        {
            ++violations;
            return false;
        }
        gh G = grid(g);
        f.kind = 1;
        f.g = g;
        f.si = Ic-G.oi;
        f.sj = Jc-G.oj;
        f.ox = (I-2*Ic==0) ? -1 : 1;
        f.oy = (J-2*Jc==0) ? -1 : 1;
        return true;
    };

    for(int l=1; l<=maxlev; ++l)
    {
        // ---- ghost cells
        vector<vector<int>> req(np), srv;
        vector<vector<int>> dst(np);

        for(int n : lev[l])
        {
            sflow_amr_patch *c = P[n];
            lexer *pp = c->pp;
            c->fill.clear();

            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                if(ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny)
                continue;

                if(pp->flagslice4[lij(pp,ii,jj)]<0)
                continue;

                int I = ii-EXT+c->I0, J = jj-EXT+c->J0;
                int o = owner(fsh(I,l),fsh(J,l));
                if(o<0)
                continue;

                sflow_amr_fill f;
                f.di=ii; f.dj=jj; f.g=-1; f.si=f.sj=0; f.ox=f.oy=0; f.slot=-1;
                f.depth = c->b->depth(ii,jj);

                if(o==me)
                {
                    if(source(l,I,J,f))
                    c->fill.push_back(f);
                    continue;
                }

                f.kind = 2;
                c->fill.push_back(f);
                req[o].push_back(I); req[o].push_back(J);
                dst[o].push_back(n); dst[o].push_back((int)c->fill.size()-1);
            }
        }

        xsetup(req,srv,2);

        sflow_amr_xplan &X = gplan[l];
        X = sflow_amr_xplan();
        gserve[l].clear();
        grecv[l].clear();

        for(int r=0; r<np; ++r)
        if(!srv[r].empty())
        {
            X.speer.push_back(r);
            vector<int> items;
            for(size_t k=0; k<srv[r].size(); k+=2)
            {
                int I=srv[r][k], J=srv[r][k+1];
                sflow_amr_fill f;
                f.di=f.dj=0; f.g=-1; f.si=f.sj=0; f.ox=f.oy=0; f.slot=-1;
                f.depth = p0->wd - bed_at(l,I,J);
                f.kind = 0;
                if(!source(l,I,J,f))
                {
                    // no data: served as dry
                    f.kind = 3;
                }
                items.push_back((int)gserve[l].size());
                gserve[l].push_back(f);
            }
            X.sitem.push_back(items);
            X.sbuf.push_back(vector<double>());
        }

        // received cells: peer index in si, position in slot
        for(int r=0; r<np; ++r)
        if(!req[r].empty())
        {
            int kp = (int)X.rpeer.size();
            X.rpeer.push_back(r);
            X.rcount.push_back((int)req[r].size()/2);
            for(size_t k=0; k<dst[r].size(); k+=2)
            {
                sflow_amr_fill *f = &P[dst[r][k]]->fill[dst[r][k+1]];
                f->si = kp;
                f->slot = (int)(k/2);
                grecv[l].push_back(f);
            }
        }

        // ---- restriction targets
        for(int n : lev[l])
        {
            sflow_amr_patch *c = P[n];
            int nb = (c->nx/2)*(c->ny/2);
            c->rgrid.assign(nb,-2); c->ric.assign(nb,0); c->rjc.assign(nb,0);

            for(int bi=0; bi<c->nx/2; ++bi)
            for(int bj=0; bj<c->ny/2; ++bj)
            {
                int k = bi*(c->ny/2)+bj;
                int Ic = (c->I0>>1)+bi, Jc = (c->J0>>1)+bj;
                if(flag0(fsh(Ic,l-1),fsh(Jc,l-1))<0)
                continue;
                int g = patch_at(l-1,Ic,Jc);
                if(g<-1)
                {
                    ++violations;
                    continue;
                }
                gh G = grid(g);
                c->rgrid[k]=g;
                c->ric[k]=Ic-G.oi;
                c->rjc[k]=Jc-G.oj;
            }
        }

        // ---- flux matching: coarse faces next to the patch boundary
        vector<vector<int>> freq(np), fsrv;
        fsend[l].clear();
        vector<vector<int>> fitems(np);

        for(int n : lev[l])
        {
            sflow_amr_patch *c = P[n];

            for(int side=0; side<4; ++side)
            {
                int ns = (side<2) ? c->ny/2 : c->nx/2;
                for(int k=0; k<ns; ++k)
                {
                    int Ic,Jc,Ii,Ji;     // outside coarse cell, inside coarse cell
                    if(side==0) { Ic=(c->I0>>1)-1; Jc=(c->J0>>1)+k; Ii=Ic+1; Ji=Jc; }
                    if(side==1) { Ic=(c->I1>>1)+1; Jc=(c->J0>>1)+k; Ii=Ic-1; Ji=Jc; }
                    if(side==2) { Ic=(c->I0>>1)+k; Jc=(c->J0>>1)-1; Ii=Ic; Ji=Jc+1; }
                    if(side==3) { Ic=(c->I0>>1)+k; Jc=(c->J1>>1)+1; Ii=Ic; Ji=Jc-1; }

                    if(flag0(fsh(Ic,l-1),fsh(Jc,l-1))<0 || flag0(fsh(Ii,l-1),fsh(Ji,l-1))<0)
                    continue;

                    int o = owner(fsh(Ic,l-1),fsh(Jc,l-1));
                    if(o<0)
                    continue;

                    if(o==me)
                    {
                        if(patch_at(l,2*Ic,2*Jc)>=0)
                        continue;

                        int g = patch_at(l-1,Ic,Jc);
                        if(g<-1)
                        {
                            ++violations;
                            continue;
                        }
                        gh G = grid(g);
                        sflow_amr_match mt;
                        mt.dir = side<2 ? 0 : 1;
                        mt.fi = Ic-G.oi - (side==1 ? 1 : 0);
                        mt.fj = Jc-G.oj - (side==3 ? 1 : 0);
                        mt.child = n;
                        mt.side = side;
                        mt.r = 2*k;
                        for(int q=0;q<6;++q) mt.val[q]=0.0;
                        match[g+1].push_back(mt);
                        continue;
                    }

                    freq[o].push_back(side); freq[o].push_back(Ic); freq[o].push_back(Jc);
                    fitems[o].push_back(n); fitems[o].push_back(side); fitems[o].push_back(2*k);
                }
            }
        }

        xsetup(freq,fsrv,3);

        sflow_amr_xplan &F = fplan[l];
        F = sflow_amr_xplan();
        frecv[l].clear();

        // I send my fine faces to the ranks of the outside coarse cells
        for(int r=0; r<np; ++r)
        if(!freq[r].empty())
        {
            F.speer.push_back(r);
            vector<int> items;
            for(size_t k=0; k<fitems[r].size(); k+=3)
            {
                items.push_back((int)fsend[l].size()/3);
                fsend[l].push_back(fitems[r][k]);
                fsend[l].push_back(fitems[r][k+1]);
                fsend[l].push_back(fitems[r][k+2]);
            }
            F.sitem.push_back(items);
            F.sbuf.push_back(vector<double>());
        }

        // and receive the fine faces of other ranks next to my coarse cells
        for(int r=0; r<np; ++r)
        if(!fsrv[r].empty())
        {
            F.rpeer.push_back(r);
            F.rcount.push_back((int)fsrv[r].size()/3);

            for(size_t k=0; k<fsrv[r].size(); k+=3)
            {
                int side=fsrv[r][k], Ic=fsrv[r][k+1], Jc=fsrv[r][k+2];
                int tgt=-1, idx=-1;

                if(patch_at(l,2*Ic,2*Jc)<0)
                {
                    int g = patch_at(l-1,Ic,Jc);
                    if(g>=-1)
                    {
                        gh G = grid(g);
                        sflow_amr_match mt;
                        mt.dir = side<2 ? 0 : 1;
                        mt.fi = Ic-G.oi - (side==1 ? 1 : 0);
                        mt.fj = Jc-G.oj - (side==3 ? 1 : 0);
                        mt.child = -1;
                        mt.side = side;
                        mt.r = 0;
                        for(int q=0;q<6;++q) mt.val[q]=0.0;
                        tgt = g+1;
                        idx = (int)rmatch[g+1].size();
                        rmatch[g+1].push_back(mt);
                    }
                    else
                    ++violations;
                }
                frecv[l].push_back(tgt);
                frecv[l].push_back(idx);
            }
        }
    }

    violations = pgc->globalisum(violations);
    if(violations>0 && p0->mpirank==0)
    cout<<"SFLOW AMR: "<<violations<<" cells without a coarser grid (nesting)"<<endl;
}

// --------------------------------------------------------------------- MPI helpers
// req[r]: keys (K ints per item) this rank wants from rank r; srv[r]: keys rank r wants from this rank
void sflow_amr::xsetup(vector<vector<int>> &req, vector<vector<int>> &srv, int K)
{
    const int np = p0->mpi_size;
    srv.assign(np,vector<int>());

    if(np==1)
    return;

    vector<int> sc(np), rc(np), sd(np,0), rd(np,0);
    for(int r=0; r<np; ++r)
    sc[r] = (int)req[r].size();

    MPI_Alltoall(&sc[0],1,MPI_INT,&rc[0],1,MPI_INT,MPI_COMM_WORLD);

    int ns=0, nr=0;
    for(int r=0; r<np; ++r)
    {
        sd[r]=ns; ns+=sc[r];
        rd[r]=nr; nr+=rc[r];
    }

    vector<int> sb(MAX(ns,1)), rb(MAX(nr,1));
    for(int r=0; r<np; ++r)
    for(int k=0; k<sc[r]; ++k)
    sb[sd[r]+k] = req[r][k];

    MPI_Alltoallv(&sb[0],&sc[0],&sd[0],MPI_INT,&rb[0],&rc[0],&rd[0],MPI_INT,MPI_COMM_WORLD);

    for(int r=0; r<np; ++r)
    srv[r].assign(rb.begin()+rd[r],rb.begin()+rd[r]+rc[r]);
}

void sflow_amr::xrun(sflow_amr_xplan &X, int nv, int tag)
{
    vector<MPI_Request> rq;
    X.rbuf.resize(X.rpeer.size());

    for(size_t k=0; k<X.rpeer.size(); ++k)
    {
        X.rbuf[k].resize((size_t)X.rcount[k]*nv);
        if(!X.rbuf[k].empty())
        {
            rq.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(&X.rbuf[k][0],(int)X.rbuf[k].size(),MPI_DOUBLE,X.rpeer[k],tag,MPI_COMM_WORLD,&rq.back());
        }
    }

    for(size_t k=0; k<X.speer.size(); ++k)
    if(!X.sbuf[k].empty())
    {
        rq.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&X.sbuf[k][0],(int)X.sbuf[k].size(),MPI_DOUBLE,X.speer[k],tag,MPI_COMM_WORLD,&rq.back());
    }

    if(!rq.empty())
    MPI_Waitall((int)rq.size(),&rq[0],MPI_STATUSES_IGNORE);
}

// --------------------------------------------------------------------- stage data
void sflow_amr::cache_stage(int s)
{
    sWL.assign(P.size()+1,nullptr);
    sUH.assign(P.size()+1,nullptr);
    sVH.assign(P.size()+1,nullptr);
    sWH.assign(P.size()+1,nullptr);

    slice *WLo,*UHo,*VHo,*WHo;
    pmom0->stage_io(s,b0,sWL[0],sUH[0],sVH[0],WLo,UHo,VHo);
    pmom0->stage_io_w(s,b0,sWH[0],WHo);

    for(size_t n=0; n<P.size(); ++n)
    {
    P[n]->pmom->stage_io(s,P[n]->b,sWL[n+1],sUH[n+1],sVH[n+1],WLo,UHo,VHo);
    P[n]->pmom->stage_io_w(s,P[n]->b,sWH[n+1],WHo);
    }
}

// well balanced prolongation from grid g (coarse cell ic,jc): the surface eta is interpolated
// with minmod slopes (switched off next to dry or solid cells), the fine depth follows from the
// fine bed; momentum through u.  Thin films (h_c < 3 A 244) copy the parent depth.
// v: WL, UH, VH, eta, U, V, wet, deep, WH, W
void sflow_amr::prolong(int g, int ic, int jc, int ox, int oy, double depthf, double *v)
{
    gh G = grid(g);
    fdm2D *pb = G.b;
    lexer *q = G.q;
    slice &WLp = *sWL[g+1], &UHp = *sUH[g+1], &VHp = *sVH[g+1], &WHp = *sWH[g+1];
    const double wd = p0->A244;

    auto get = [&](int a, int bb, double &e, double &u, double &vv, double &ww, int &wt)
    {
        double wlc = WLp(a,bb);
        e = wlc - pb->depth(a,bb);
        wt = q->wet[lij(q,a,bb)];
        if(q->flagslice4[lij(q,a,bb)]<0)
        wt = -1;
        double wlvl = fabs(wlc)>wd ? wlc : 1.0e20;
        u = wt==1 ? UHp(a,bb)/wlvl : 0.0;
        vv = wt==1 ? VHp(a,bb)/wlvl : 0.0;
        ww = wt==1 ? WHp(a,bb)/wlvl : 0.0;
    };

    auto dry = [&]()
    {
        v[0]=wd; v[1]=v[2]=0.0; v[3]=wd-depthf; v[4]=v[5]=0.0; v[6]=0.0; v[7]=0.0; v[8]=v[9]=0.0;
    };

    double e0,u0,v0,w0v; int w0;
    get(ic,jc,e0,u0,v0,w0v,w0);

    if(w0!=1)
    {
        dry();
        return;
    }

    auto fin = [&](double wl, double uh, double vh, double wh)
    {
        double wlvl = fabs(wl)>wd ? wl : 1.0e20;
        v[0]=wl; v[1]=uh; v[2]=vh;
        v[3]=wl-depthf;
        v[4]=uh/wlvl;
        v[5]=vh/wlvl*p0->y_dir;
        v[6]=1.0; v[7]=1.0;
        v[8]=wh; v[9]=wh/wlvl;
    };

    double hc = WLp(ic,jc);
    if(hc<3.0*wd)
    {
        fin(hc,UHp(ic,jc),VHp(ic,jc),WHp(ic,jc));
        return;
    }

    double eE,uE,vE,wEv,eW,uW,vW,wWv,eN,uN,vN,wNv,eS,uS,vS,wSv;
    int wE,wW,wN,wS;
    get(ic+1,jc,eE,uE,vE,wEv,wE);
    get(ic-1,jc,eW,uW,vW,wWv,wW);
    get(ic,jc+1,eN,uN,vN,wNv,wN);
    get(ic,jc-1,eS,uS,vS,wSv,wS);

    double sxe=0.0,sye=0.0,sxu=0.0,syu=0.0,sxv=0.0,syv=0.0,sxw=0.0,syw=0.0;
    if(wE==1 && wW==1 && wN==1 && wS==1)
    {
        sxe = mmod(eE-e0,e0-eW); sye = mmod(eN-e0,e0-eS);
        sxu = mmod(uE-u0,u0-uW); syu = mmod(uN-u0,u0-uS);
        sxv = mmod(vE-v0,v0-vW); syv = mmod(vN-v0,v0-vS);
        sxw = mmod(wEv-w0v,w0v-wWv); syw = mmod(wNv-w0v,w0v-wSv);
    }

    const double fx = 0.25*ox, fy = 0.25*oy;
    double ef = e0 + sxe*fx + sye*fy;
    double hf = ef + depthf;

    if(hf<=wd+eps)
    {
        dry();
        return;
    }

    fin(hf, hf*(u0 + sxu*fx + syu*fy), hf*(v0 + sxv*fx + syv*fy), hf*(w0v + sxw*fx + syw*fy));
}

void sflow_amr::eval_fill(const sflow_amr_fill &f, double *v)
{
    if(f.kind==1)
    {
        prolong(f.g,f.si,f.sj,f.ox,f.oy,f.depth,v);
        return;
    }

    if(f.kind==0)
    {
        gh G = grid(f.g);
        v[0] = (*sWL[f.g+1])(f.si,f.sj);
        v[1] = (*sUH[f.g+1])(f.si,f.sj);
        v[2] = (*sVH[f.g+1])(f.si,f.sj);
        v[3] = G.b->eta(f.si,f.sj);
        v[4] = G.b->U(f.si,f.sj);
        v[5] = G.b->V(f.si,f.sj);
        v[6] = G.q->wet[lij(G.q,f.si,f.sj)];
        v[7] = G.q->deep[lij(G.q,f.si,f.sj)];
        v[8] = (*sWH[f.g+1])(f.si,f.sj);
        v[9] = G.b->W(f.si,f.sj);
        return;
    }

    const double wd = p0->A244;
    v[0]=wd; v[1]=v[2]=0.0; v[3]=wd-f.depth; v[4]=v[5]=0.0; v[6]=v[7]=0.0; v[8]=v[9]=0.0;
}

// cells of the level-l patches outside their interior, then the wall conditions
void sflow_amr::fill_level(ghostcell *pgc, int l, int s)
{
    sflow_amr_xplan &X = gplan[l];
    double v[NV];

    // cells served to other ranks
    for(size_t k=0; k<X.speer.size(); ++k)
    {
        vector<double> &sb = X.sbuf[k];
        sb.resize(X.sitem[k].size()*NV);
        for(size_t m=0; m<X.sitem[k].size(); ++m)
        eval_fill(gserve[l][X.sitem[k][m]],&sb[m*NV]);
    }

    xrun(X,NV,7100+l);

    auto store = [&](sflow_amr_patch *c, int n, int ii, int jj, const double *w)
    {
        lexer *pp = c->pp;
        fdm2D *b = c->b;
        (*sWL[n+1])(ii,jj) = w[0];
        (*sUH[n+1])(ii,jj) = w[1];
        (*sVH[n+1])(ii,jj) = w[2];
        b->eta(ii,jj) = w[3];
        b->U(ii,jj) = w[4];
        b->V(ii,jj) = w[5];
        b->hp(ii,jj) = w[0];
        pp->wet[lij(pp,ii,jj)] = (int)w[6];
        pp->deep[lij(pp,ii,jj)] = (int)w[7];
        (*sWH[n+1])(ii,jj) = w[8];
        b->W(ii,jj) = w[9];
    };

    comms_off guard(pgc);

    for(int n : lev[l])
    {
        sflow_amr_patch *c = P[n];
        for(auto &f : c->fill)
        {
            if(f.kind==2)
            {
                int k = f.si, m = f.slot;
                store(c,n,f.di,f.dj,&X.rbuf[k][(size_t)m*NV]);
                continue;
            }
            eval_fill(f,v);
            store(c,n,f.di,f.dj,v);
        }
    }

    for(int n : lev[l])
    apply_bc(pgc,*P[n],s);
}

// wall ghost cells of a patch from its (filled) fluid cells, as at the end of a level-0 stage
void sflow_amr::apply_bc(ghostcell *pgc, sflow_amr_patch &c, int s)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;
    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);

    slice *WHi,*WHo;
    c.pmom->stage_io_w(s,b,WHi,WHo);

    int gcval_eta = 50+p0->F50;

    pgc->gcsl_start4(pp,b->eta,gcval_eta);
    pgc->gcsl_start4Vint(pp,pp->wet,50);
    pgc->gcsl_start4Vint(pp,pp->deep,50);
    c.pmom->ghostcells(pp,b,pgc,*UHi,*VHi,*WHi,*WLi);
}

void sflow_amr::restrict_patch(lexer *p, sflow_amr_patch &c, int s)
{
    fdm2D *b = c.b;
    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);
    slice *WHi,*WHo;
    c.pmom->stage_io_w(s,b,WHi,WHo);

    const int nby = c.ny/2;

    for(int bi=0; bi<c.nx/2; ++bi)
    for(int bj=0; bj<nby; ++bj)
    {
        int k = bi*nby+bj;
        int g = c.rgrid[k];
        if(g<-1)
        continue;

        gh G = grid(g);
        slice *pWLi,*pUHi,*pVHi,*pWLo,*pUHo,*pVHo;
        G.m->stage_io(s,G.b,pWLi,pUHi,pVHi,pWLo,pUHo,pVHo);
        slice *pWHi,*pWHo;
        G.m->stage_io_w(s,G.b,pWHi,pWHo);

        int ic = c.ric[k], jc = c.rjc[k];
        int i0 = EXT+2*bi, j0 = EXT+2*bj;

        double wl = 0.25*((*WLo)(i0,j0)+(*WLo)(i0+1,j0)+(*WLo)(i0,j0+1)+(*WLo)(i0+1,j0+1));
        double uh = 0.25*((*UHo)(i0,j0)+(*UHo)(i0+1,j0)+(*UHo)(i0,j0+1)+(*UHo)(i0+1,j0+1));
        double vh = 0.25*((*VHo)(i0,j0)+(*VHo)(i0+1,j0)+(*VHo)(i0,j0+1)+(*VHo)(i0+1,j0+1));
        double wh = 0.25*((*WHo)(i0,j0)+(*WHo)(i0+1,j0)+(*WHo)(i0,j0+1)+(*WHo)(i0+1,j0+1));

        (*pWLo)(ic,jc)=wl; (*pUHo)(ic,jc)=uh; (*pVHo)(ic,jc)=vh; (*pWHo)(ic,jc)=wh;

        int w = wl>p->A244+eps ? 1 : 0;
        G.q->wet[lij(G.q,ic,jc)] = w;

        double wlvl = fabs(wl)>p->A244 ? wl : 1.0e20;
        G.b->eta(ic,jc) = wl - G.b->depth(ic,jc);
        G.b->U(ic,jc) = w==1 ? uh/wlvl : 0.0;
        G.b->V(ic,jc) = w==1 ? vh/wlvl : 0.0;
        G.b->W(ic,jc) = w==1 ? wh/wlvl : 0.0;
        G.b->hp(ic,jc) = wl;
    }
}

// level-0 halo after the restriction: the neighbour ranks see the restricted values
void sflow_amr::exchange_level0(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(p->mpi_size<=1)
    return;

    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    pmom0->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);

    slice *WHi,*WHo;
    pmom0->stage_io_w(s,b,WHi,WHo);

    slice *f[9] = {WLo,UHo,VHo,WHo,&b->eta,&b->U,&b->V,&b->W,&b->hp};
    for(int k=0; k<9; ++k)
    {
        pgc->gcslparax(p,*f[k],4);
        pgc->gcslparacox(p,*f[k],10);
    }

    pgc->gcslparaxV_int(p,p->wet,4);
    pgc->gcslparacoxV_int(p,p->wet,50);
}

// fine face values of the level-l patches on the partition edges, to the rank of the coarse cell
void sflow_amr::exchange_fluxes(int l)
{
    sflow_amr_xplan &F = fplan[l];

    for(size_t k=0; k<F.speer.size(); ++k)
    {
        vector<double> &sb = F.sbuf[k];
        sb.resize(F.sitem[k].size()*5);
        for(size_t m=0; m<F.sitem[k].size(); ++m)
        {
            int it = F.sitem[k][m];
            sflow_amr_patch *c = P[fsend[l][3*it]];
            int side = fsend[l][3*it+1], r = fsend[l][3*it+2];
            for(int ip=0; ip<5; ++ip)
            sb[m*5+ip] = 0.5*(c->rec[ip][side][r]+c->rec[ip][side][r+1]);
        }
    }

    xrun(F,5,7200+l);

    size_t it=0;
    for(size_t k=0; k<F.rpeer.size(); ++k)
    for(int m=0; m<F.rcount[k]; ++m)
    {
        int tgt = frecv[l][2*it], idx = frecv[l][2*it+1];
        if(tgt>=0)
        for(int ip=0; ip<5; ++ip)
        rmatch[tgt][idx].val[ip] = F.rbuf[k][(size_t)m*5+ip];
        ++it;
    }
}

// --------------------------------------------------------------------- stages
void sflow_amr::step_begin(lexer *p, fdm2D *b, ghostcell *pgc)
{
    comms_off guard(pgc);

    for(auto c : P)
    {
        c->pp->dt = p->dt;
        c->pp->dt_old = p->dt_old;
        c->pp->simtime = p->simtime;
        c->pp->count = p->count;
        c->pmom->inflow(c->pp,c->b,pgc,pflow_void);
    }

    // body on the patches: X 10 3 updates the ship pressure once per step (as level 0),
    // X 10 2 moves the body in every stage
    ship_patches(shipmode==3);
}

void sflow_amr::stage_begin(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(maxlev<1 || patches_total==0)
    return;

    double t0 = MPI_Wtime();
    cache_stage(s);

    for(int l=1; l<=maxlev; ++l)
    fill_level(pgc,l,s);

    if(shipmode==2 && s>0)
    ship_patches(false);
    tm[0] += MPI_Wtime()-t0;


    // finest first: the coarser grids take the recorded fine fluxes
    for(int l=maxlev; l>=1; --l)
    {
        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        P[n]->pmom->rk_stage(P[n]->pp,P[n]->b,pgc,s);
        }
        tm[1] += MPI_Wtime()-t0;

        t0 = MPI_Wtime();
        exchange_fluxes(l);
        tm[2] += MPI_Wtime()-t0;
    }
}

void sflow_amr::stage_end(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(maxlev<1 || patches_total==0)
    return;

    double t0 = MPI_Wtime();
    {
    comms_off guard(pgc);
    for(int l=maxlev; l>=1; --l)
    for(int n : lev[l])
    restrict_patch(p,*P[n],s);
    }

    exchange_level0(p,b,pgc,s);
    tm[3] += MPI_Wtime()-t0;
}

void sflow_amr::step_end(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // end of the step as for level 0 in sflow_f::mainloop: breaking flags cleared (A 248 0),
    // water depth and wet-dry state updated
    {
    comms_off guard(pgc);
    for(auto c : P)
    {
        lexer *pp = c->pp;
        c->pmom->rk_finish(pp,c->b,pgc);

        if(p->A248==0)
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        c->b->breaking(ii,jj)=0;

        c->pfsf->depth_update(pp,c->b,pgc,c->b->WL);
    }
    }

    if(regrid_int>0 && p->count%regrid_int==0)
    {
    double t0 = MPI_Wtime();
    regrid(p,b,pgc,false);
    tm[4] += MPI_Wtime()-t0;
    }
}

// --------------------------------------------------------------------- flux matching
void sflow_amr::hll_hook(lexer *p, fdm2D *b, int ipol, int id)
{
    slice &Fx = (ipol==4) ? b->FEx : b->Fx;
    slice &Fy = (ipol==4) ? b->FEy : b->Fy;

    // a patch records its boundary faces
    if(id>0)
    {
        sflow_amr_patch &c = *P[id-1];
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EXT-1, jh = EXT+c.ny-1;

        for(int r=0; r<c.ny; ++r)
        {
            c.rec[ipol][0][r] = Fx(il,EXT+r);
            c.rec[ipol][1][r] = Fx(ih,EXT+r);
        }
        for(int r=0; r<c.nx; ++r)
        {
            c.rec[ipol][2][r] = Fy(EXT+r,jl);
            c.rec[ipol][3][r] = Fy(EXT+r,jh);
        }

        if(ipol==4)
        {
            for(int r=0; r<c.ny; ++r)
            {
                c.rec[0][0][r] = b->dfx(il,EXT+r);
                c.rec[0][1][r] = b->dfx(ih,EXT+r);
            }
            for(int r=0; r<c.nx; ++r)
            {
                c.rec[0][2][r] = b->dfy(EXT+r,jl);
                c.rec[0][3][r] = b->dfy(EXT+r,jh);
            }
        }
    }

    if(id>=(int)match.size())
    return;

    // coarse faces next to the patches of this rank: mean of the two fine faces
    for(auto &m : match[id])
    {
        sflow_amr_patch &c = *P[m.child];
        double val = 0.5*(c.rec[ipol][m.side][m.r]+c.rec[ipol][m.side][m.r+1]);

        if(m.dir==0)
        Fx(m.fi,m.fj) = val;
        else
        Fy(m.fi,m.fj) = val;

        if(ipol==4)
        {
            double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+1]);
            if(m.dir==0)
            b->dfx(m.fi,m.fj) = df;
            else
            b->dfy(m.fi,m.fj) = df;
        }
    }

    // coarse faces next to patches of other ranks
    for(auto &m : rmatch[id])
    {
        if(m.dir==0)
        Fx(m.fi,m.fj) = m.val[ipol];
        else
        Fy(m.fi,m.fj) = m.val[ipol];

        if(ipol==4)
        {
            if(m.dir==0)
            b->dfx(m.fi,m.fj) = m.val[0];
            else
            b->dfy(m.fi,m.fj) = m.val[0];
        }
    }
}

// --------------------------------------------------------------------- time step
void sflow_amr::timestep(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // initial time step (sflow_etimestep::ini): linear in DXM, so scale with the finest patch
    if(p->count==0)
    {
        double r=1.0;
        for(auto c : P)
        r = MIN(r, c->pp->DXM/p->DXM);

        r = pgc->globalmin(r);
        p->dt *= r;
        p->dt_old = p->dt;
        return;
    }

    // same CFL as sflow_etimestep, over the wet real patch cells
    const double g = fabs(p->W22);
    double cmin = 1.0e20;

    for(auto c : P)
    {
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(pp->wet[lij(pp,ii,jj)]==1 && pp->flagslice4[lij(pp,ii,jj)]>0)
        {
            double cc = sqrt(g*MAX(pb->WL(ii,jj),p->A244));
            cmin = MIN(cmin, pp->DXN[ii+marge]/(fabs(pb->U(ii,jj))+cc));
            cmin = MIN(cmin, pp->DYN[jj+marge]/(fabs(pb->V(ii,jj))+cc));

            if(p->A219==2)
            cmin = MIN(cmin, pp->DXN[ii+marge]/(fabs(pb->U(ii,jj))>1.0e-20?fabs(pb->U(ii,jj)):1.0e-20));
        }
    }

    double dtp = p->N47*2.0*cmin;
    dtp = pgc->globalmin(dtp);

    if(p->N48==1)
    p->dt = MIN(p->dt,dtp);
}

// --------------------------------------------------------------------- output
double sflow_amr::mass(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double m=0.0;
    SLICELOOP4
    m += b->WL(i,j)*p->DXN[IP]*p->DYN[JP];

    return pgc->globalsum(m);
}

void sflow_amr::print(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    bool doprint=false;

    if((p->count%p->P181==0 && p->P182<0.0 && p->P10==1 && p->P181>0) || (p->count==0 && p->P182<0.0))
    doprint=true;

    if((p->simtime>printtime_amr && p->P182>0.0 && p->P10==1) || (p->count==0 && p->P182>0.0))
    {
        doprint=true;
        printtime_amr += p->P182;
    }

    double m = mass(p,b,pgc);

    if(p->P51>0)
    gauges(p,b,pgc);

    if(p->mpirank==0 && (p->count%p->P12==0 || doprint))
    logout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<p->dt<<" \t "<<patches_total<<" \t "<<cells_total<<" \t "
          <<setprecision(15)<<m<<" \t "<<setprecision(6)<<(m-m0)/(fabs(m0)>1.0e-20?m0:1.0)<<endl;

    if(p->mpirank==0 && doprint)
    cout<<"SFLOW AMR: "<<patches_total<<" patches, "<<cells_total<<" cells; time in patch stages "<<setprecision(4)<<tm[1]<<" s, fill "<<tm[0]<<" s, restriction "<<tm[3]<<" s, regrid "<<tm[4]<<" s";
    if(p->mpirank==0 && doprint && shipmode>0)
    cout<<", body "<<tm[8]<<" s";
    if(p->mpirank==0 && doprint && nh==1)
    cout<<"; pressure "<<tm[5]<<" s (preconditioner "<<tm[6]<<" s, operator "<<tm[7]<<" s), mean iterations "<<(nh_solves>0 ? double(nh_it_total)/nh_solves : 0.0);
    if(p->mpirank==0 && doprint)
    cout<<endl;

    if(!doprint)
    return;

    write_vtr0(p,b);

    for(int n=0; n<(int)P.size(); ++n)
    write_vtr(p,*P[n],n);

    // multiblock index: level 0 of every rank + all patches
    int np = (int)P.size();
    vector<int> all(p->mpi_size,0);
    MPI_Allgather(&np,1,MPI_INT,&all[0],1,MPI_INT,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        char name[256];
        snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i.vtm",printcount_amr);
        ofstream out(name);
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n<vtkMultiBlockDataSet>\n";
        out<<"<Block index=\"0\" name=\"level 0\">\n";
        for(int r=0; r<p->mpi_size; ++r)
        out<<"<DataSet index=\""<<r<<"\" file=\"REEF3D-SFLOW-AMR-L0-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<".vtr\"/>\n";
        out<<setfill(' ')<<"</Block>\n<Block index=\"1\" name=\"patches\">\n";
        int idx=0;
        for(int r=0; r<p->mpi_size; ++r)
        for(int q=0; q<all[r]; ++q)
        {
        out<<"<DataSet index=\""<<idx<<"\" file=\"REEF3D-SFLOW-AMR-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<"-"<<setw(4)<<q+1<<".vtr\"/>\n";
        out<<setfill(' ');
        ++idx;
        }
        out<<"</Block>\n</vtkMultiBlockDataSet>\n</VTKFile>\n";
        out.close();
    }

    ++printcount_amr;
}

// surface elevation at the P 51 gauges from the finest grid that holds them
void sflow_amr::gauges(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const int ng = p->P51;
    vector<double> lv(ng,-1.0), val(ng,0.0);

    auto find = [&](lexer *q, int i0, int i1, int j0, int j1, double x, double y, int &ii, int &jj)
    {
        ii=-1; jj=-1;
        for(int a=i0; a<=i1; ++a)
        if(x>=q->XN[a+marge] && x<q->XN[a+1+marge]) { ii=a; break; }
        for(int a=j0; a<=j1; ++a)
        if(y>=q->YN[a+marge] && y<q->YN[a+1+marge]) { jj=a; break; }
        return ii>=0 && jj>=0;
    };

    for(int k=0; k<ng; ++k)
    {
        int ii,jj;
        if(find(p,0,NX0-1,0,NY0-1,p->P51_x[k],p->P51_y[k],ii,jj))
        {
            lv[k]=0.0;
            val[k]=b->eta(ii,jj);
        }
        for(auto c : P)
        if(c->lev>lv[k] && find(c->pp,EXT,EXT+c->nx-1,EXT,EXT+c->ny-1,p->P51_x[k],p->P51_y[k],ii,jj))
        {
            lv[k]=c->lev;
            val[k]=c->b->eta(ii,jj);
        }
    }

    vector<double> lmax(ng);
    MPI_Allreduce(&lv[0],&lmax[0],ng,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    for(int k=0; k<ng; ++k)
    if(lv[k]!=lmax[k])
    val[k]=0.0;
    vector<double> vs(ng);
    MPI_Allreduce(&val[0],&vs[0],ng,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        if(!gaugeout.is_open())
        {
            gaugeout.open("./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_gauges.dat");
            gaugeout<<"# simtime";
            for(int k=0; k<ng; ++k)
            gaugeout<<" \t eta("<<p->P51_x[k]<<","<<p->P51_y[k]<<")";
            gaugeout<<" \t levels"<<endl;
        }
        gaugeout<<setprecision(10)<<p->simtime;
        for(int k=0; k<ng; ++k)
        gaugeout<<" \t "<<setprecision(10)<<vs[k];
        for(int k=0; k<ng; ++k)
        gaugeout<<" \t "<<int(lmax[k]);
        gaugeout<<endl;
    }
}

// level 0 of this rank, cell centred, same layout as the patches
void sflow_amr::write_vtr0(lexer *p, fdm2D *b)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-L0-%08i-%04i.vtr",printcount_amr,p->mpirank+1);

    const int nx=p->knox, ny=p->knoy, m=marge;
    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<RectilinearGrid WholeExtent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">0</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=0; jj<ny; ++jj)
        {
            for(int ii=0; ii<nx; ++ii)
            out<<setprecision(17)<<f(ii,jj)+shift<<" ";
            out<<"\n";
        }
        out<<"</DataArray>\n";
    };
    field("eta",b->eta,0.0);
    field("elevation",b->eta,p->wd);
    field("WL",b->WL,0.0);
    field("u",b->U,0.0);
    field("v",b->V,0.0);
    field("bed",b->bed,0.0);
    field("press",b->press,0.0);
    field("w",b->W,0.0);

    out<<"<DataArray type=\"Int32\" Name=\"wetdry\" format=\"ascii\">\n";
    for(int jj=0; jj<ny; ++jj)
    {
        for(int ii=0; ii<nx; ++ii)
        out<<p->wet[lij(p,ii,jj)]<<" ";
        out<<"\n";
    }
    out<<"</DataArray>\n</CellData>\n<Coordinates>\n";
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=0; ii<=nx; ++ii) out<<setprecision(12)<<p->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=0; jj<=ny; ++jj) out<<setprecision(12)<<p->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}

void sflow_amr::write_vtr(lexer *p, sflow_amr_patch &c, int id)
{
    lexer *pp = c.pp;
    fdm2D *pb = c.b;
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i-%04i-%04i.vtr",printcount_amr,p->mpirank+1,id+1);

    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<RectilinearGrid WholeExtent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">"<<c.lev<<"</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=EXT; jj<EXT+c.ny; ++jj)
        {
            for(int ii=EXT; ii<EXT+c.nx; ++ii)
            out<<setprecision(17)<<f(ii,jj)+shift<<" ";
            out<<"\n";
        }
        out<<"</DataArray>\n";
    };
    field("eta",pb->eta,0.0);
    field("elevation",pb->eta,p->wd);
    field("WL",pb->WL,0.0);
    field("u",pb->U,0.0);
    field("v",pb->V,0.0);
    field("bed",pb->bed,0.0);
    field("press",pb->press,0.0);
    field("w",pb->W,0.0);

    out<<"<DataArray type=\"Int32\" Name=\"wetdry\" format=\"ascii\">\n";
    for(int jj=EXT; jj<EXT+c.ny; ++jj)
    {
        for(int ii=EXT; ii<EXT+c.nx; ++ii)
        out<<pp->wet[lij(pp,ii,jj)]<<" ";
        out<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"</CellData>\n<Coordinates>\n";

    const int m = marge;
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=EXT; ii<=EXT+c.nx; ++ii) out<<setprecision(12)<<pp->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=EXT; jj<=EXT+c.ny; ++jj) out<<setprecision(12)<<pp->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}
