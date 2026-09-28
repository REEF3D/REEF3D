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
#include"ioflow_void.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<cstdio>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
inline int fdiv2(int a) { return (a>=0) ? a/2 : -((1-a)/2); }
inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b)?a:b;
}
// index of cell (i,j) in the lexer 2D arrays (wet, flagslice4)
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

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
    pBC = ppBC;
    p6dof = pp6dof;
    phll0 = pphll;
    pmom0 = ppmom;
    maxlev = p->A270;
    patches_total = 0;
    printcount_amr = 0;
    printtime_amr = 0.0;
    m0 = 0.0;

    pflow_void = new ioflow_v(p,pgc,pBC);

    phll0->amr = this;
    phll0->amr_id = 0;
}

sflow_amr::~sflow_amr()
{
}

// --------------------------------------------------------------------- parent access
lexer* sflow_amr::plex(sflow_amr_patch &c) { return c.parent<0 ? p0 : P[c.parent].pp; }
fdm2D* sflow_amr::pfdm(sflow_amr_patch &c) { return c.parent<0 ? b0 : P[c.parent].b; }
sflow_momentum_RK3* sflow_amr::pmom(sflow_amr_patch &c) { return c.parent<0 ? pmom0 : P[c.parent].pmom; }
int sflow_amr::pwet(sflow_amr_patch &c, int ii, int jj) { lexer *q=plex(c); return q->wet[lij(q,ii,jj)]; }
void sflow_amr::set_pwet(sflow_amr_patch &c, int ii, int jj, int w) { lexer *q=plex(c); q->wet[lij(q,ii,jj)]=w; }

// --------------------------------------------------------------------- setup
void sflow_amr::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1 || p->A276==0)
    return;

    // v1 scope
    int ok=1;
    if(p->A220!=0) ok=0;
    if(p->A210==2) ok=0;
    if(p->A212!=0) ok=0;
    if(p->A260!=0) ok=0;
    if(p->S10!=0) ok=0;
    if(p->X10>=2) ok=0;
    if(p->W90!=0) ok=0;
    if(p->j_dir!=1) ok=0;

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"SFLOW AMR (A 270): only for A 220 0, A 210 3, A 212 0, A 260 0, S 10 0, X 10 0/1, W 90 0 and 2D grids -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    const int nest = 3;
    const int grow = 2*(maxlev-1);   // level-0 cells between the level-1 and the level-L region

    for(int q=0; q<p->A276; ++q)
    {
        // box in level-0 cells of this rank
        int ilo=p->knox, ihi=-1, jlo=p->knoy, jhi=-1;

        for(i=0; i<p->knox; ++i)
        if(p->XP[IP]>=p->A276_xs[q] && p->XP[IP]<=p->A276_xe[q])
        {
            ilo = MIN(ilo,i);
            ihi = MAX(ihi,i);
        }

        for(j=0; j<p->knoy; ++j)
        if(p->YP[JP]>=p->A276_ys[q] && p->YP[JP]<=p->A276_ye[q])
        {
            jlo = MIN(jlo,j);
            jhi = MAX(jhi,j);
        }

        // clip: the level-1 region stays nest cells away from the partition edges
        ilo = MAX(ilo, nest+grow);
        jlo = MAX(jlo, nest+grow);
        ihi = MIN(ihi, p->knox-1-nest-grow);
        jhi = MIN(jhi, p->knoy-1-nest-grow);

        if(ilo>ihi || jlo>jhi)
        continue;

        // walls and solids: none within nest cells of the level-1 region
        int clean=1;
        for(i=ilo-grow-nest; i<=ihi+grow+nest; ++i)
        for(j=jlo-grow-nest; j<=jhi+grow+nest; ++j)
        if(p->flagslice4[IJ]<0)
        clean=0;

        // chains must not overlap (v1: no siblings)
        for(auto &c : P)
        if(c.lev==1)
        if(!(ilo-grow-nest>c.ihi0 || ihi+grow+nest<c.ilo0 || jlo-grow-nest>c.jhi0 || jhi+grow+nest<c.jlo0))
        clean=0;

        if(clean==0)
        {
            cout<<"SFLOW AMR: rank "<<p->mpirank<<", A 276 box "<<q+1<<" skipped (wall, solid or overlap within "<<nest<<" cells)"<<endl;
            continue;
        }

        int parent=-1;
        for(int l=1; l<=maxlev; ++l)
        {
            int g = 2*(maxlev-l);
            make_patch(p,b,pgc,l,parent,ilo-g,ihi+g,jlo-g,jhi+g);
            parent = (int)P.size()-1;
        }
    }

    // coarse to fine
    order.clear();
    for(int l=1; l<=maxlev; ++l)
    for(int n=0; n<(int)P.size(); ++n)
    if(P[n].lev==l)
    order.push_back(n);

    {
    comms_off guard(pgc);

    for(int n : order)
    {
        ini_bed(p,P[n]);
        ini_state(p,pgc,P[n]);
    }
    }

    // make level 0 consistent with the patches (conservative averages)
    for(int s=2; s<3; ++s)
    {
    comms_off guard(pgc);
    for(int k=(int)order.size()-1; k>=0; --k)
    restrict_patch(p,P[order[k]],s);
    }

    patches_total = pgc->globalisum((int)P.size());
    m0 = mass(p,b,pgc);

    mkdir("./REEF3D_SFLOW_AMR",0777);
    
    if(p->mpirank==0)
    {
        logout.open("./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_log.dat");
        logout<<"# count \t simtime \t dt \t patches \t mass \t rel. mass change"<<endl;
        cout<<"SFLOW AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es)"<<endl;
    }
}

void sflow_amr::make_patch(lexer *p, fdm2D *b, ghostcell *pgc, int lev, int parent, int ilo0, int ihi0, int jlo0, int jhi0)
{
    sflow_amr_patch c;
    c.id = (int)P.size();
    c.lev = lev;
    c.parent = parent;
    c.ilo0=ilo0; c.ihi0=ihi0; c.jlo0=jlo0; c.jhi0=jhi0;

    int r = 1<<lev;
    c.nx = (ihi0-ilo0+1)*r;
    c.ny = (jhi0-jlo0+1)*r;

    if(parent<0)
    {
        c.pi0=ilo0; c.pi1=ihi0;
        c.pj0=jlo0; c.pj1=jhi0;
    }

    if(parent>=0)
    {
        sflow_amr_patch &Q = P[parent];
        int rp = 1<<(lev-1);
        c.pi0 = 1 + (ilo0-Q.ilo0)*rp;
        c.pj0 = 1 + (jlo0-Q.jlo0)*rp;
        c.pi1 = c.pi0 + (ihi0-ilo0+1)*rp - 1;
        c.pj1 = c.pj0 + (jhi0-jlo0+1)*rp - 1;
    }

    lexer *par = (parent<0) ? p : P[parent].pp;

    c.pp = new lexer(*p,1);
    build_lexer(p,c,par);

    {
    comms_off guard(pgc);
    lexer *pp = c.pp;

    c.b = new fdm2D(pp);
    c.phll = new sflow_HLL(pp,pgc,pBC);
    c.phll->amr = this;
    c.phll->amr_id = c.id+1;
    c.pss = new sflow_signal_speed(pp);

    if(p->A211<=3)
    c.precon = new sflow_reconstruct_hires(pp,pBC);
    if(p->A211>=4)
    c.precon = new sflow_reconstruct_weno(pp,pBC,1);

    c.pdiff = new sflow_diffusion_void(pp);
    c.ppress = new sflow_hydrostatic(pp,c.b,pBC);
    c.pfsf = new sflow_eta(pp,c.b,pgc,pBC);
    c.psfdf = new sflow_forcing(pp);
    c.pmom = new sflow_momentum_RK3(pp,c.b,pgc,c.phll,c.pss,c.precon,c.pdiff,c.ppress,nullptr,nullptr,pflow_void,c.pfsf,c.psfdf,p6dof);
    }

    for(int ip=0; ip<5; ++ip)
    {
        c.rec[ip][0].assign(c.ny,0.0);
        c.rec[ip][1].assign(c.ny,0.0);
        c.rec[ip][2].assign(c.nx,0.0);
        c.rec[ip][3].assign(c.nx,0.0);
    }

    P.push_back(c);

    if(parent>=0)
    P[parent].children.push_back(c.id);
}

// patch geometry: nodes are the parent's nodes and their midpoints
void sflow_amr::build_lexer(lexer *p, sflow_amr_patch &c, lexer *par)
{
    lexer *pp = c.pp;
    const int ms = pp->margin;   // ghost layers of the slices
    const int m = marge;         // offset of the coordinate arrays (IP = i+marge)

    pp->knox = c.nx+1;
    pp->knoy = c.ny+1;
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
    const int pnxa = par->knox+1+4*m;
    const int pnya = par->knoy+1+4*m;

    pp->Darray(pp->XN,nxa); pp->Darray(pp->XP,nxa); pp->Darray(pp->DXN,nxa); pp->Darray(pp->DXP,nxa);
    pp->Darray(pp->YN,nya); pp->Darray(pp->YP,nya); pp->Darray(pp->DYN,nya); pp->Darray(pp->DYP,nya);
    pp->ZN = p->ZN; pp->ZP = p->ZP; pp->DZN = p->DZN; pp->DZP = p->DZP;

    auto pxn = [&](int kn)   // parent node kn, linear extrapolation outside the parent arrays
    {
        int a = kn+m;
        if(a<0)      return par->XN[0] + a*(par->XN[1]-par->XN[0]);
        if(a>=pnxa)  return par->XN[pnxa-1] + (a-pnxa+1)*(par->XN[pnxa-1]-par->XN[pnxa-2]);
        return par->XN[a];
    };
    auto pyn = [&](int kn)
    {
        int a = kn+m;
        if(a<0)      return par->YN[0] + a*(par->YN[1]-par->YN[0]);
        if(a>=pnya)  return par->YN[pnya-1] + (a-pnya+1)*(par->YN[pnya-1]-par->YN[pnya-2]);
        return par->YN[a];
    };

    for(int a=0; a<nxa; ++a)
    {
        int r = a-m-1;          // fine node relative to the patch origin
        int k = fdiv2(r);
        pp->XN[a] = (r-2*k==0) ? pxn(c.pi0+k) : 0.5*(pxn(c.pi0+k)+pxn(c.pi0+k+1));
    }
    for(int a=0; a<nya; ++a)
    {
        int r = a-m-1;
        int k = fdiv2(r);
        pp->YN[a] = (r-2*k==0) ? pyn(c.pj0+k) : 0.5*(pyn(c.pj0+k)+pyn(c.pj0+k+1));
    }

    for(int a=0; a<nxa-1; ++a) { pp->XP[a] = 0.5*(pp->XN[a]+pp->XN[a+1]); pp->DXN[a] = pp->XN[a+1]-pp->XN[a]; }
    for(int a=0; a<nya-1; ++a) { pp->YP[a] = 0.5*(pp->YN[a]+pp->YN[a+1]); pp->DYN[a] = pp->YN[a+1]-pp->YN[a]; }
    pp->XP[nxa-1] = pp->XP[nxa-2] + pp->DXN[nxa-2]; pp->DXN[nxa-1] = pp->DXN[nxa-2];
    pp->YP[nya-1] = pp->YP[nya-2] + pp->DYN[nya-2]; pp->DYN[nya-1] = pp->DYN[nya-2];
    for(int a=0; a<nxa-1; ++a) pp->DXP[a] = pp->XP[a+1]-pp->XP[a];
    for(int a=0; a<nya-1; ++a) pp->DYP[a] = pp->YP[a+1]-pp->YP[a];
    pp->DXP[nxa-1] = pp->DXP[nxa-2];
    pp->DYP[nya-1] = pp->DYP[nya-2];

    // mean spacing as grid::gridspacing
    double s=0.0; int n=0;
    for(int ii=0; ii<pp->knox; ++ii) { s += pp->DXP[ii+m]; ++n; }
    for(int jj=0; jj<pp->knoy; ++jj) { s += pp->DYP[jj+m]; ++n; }
    pp->DXM = s/double(n);
    pp->DXD = p->DXD;
    pp->DYD = p->DYD;
    pp->DYM = pp->DXM;

    // flags: all fluid (patches keep away from walls and solids)
    const int nsl = pp->imax*pp->jmax;
    pp->Iarray(pp->flagslice1,nsl);
    pp->Iarray(pp->flagslice2,nsl);
    pp->Iarray(pp->flagslice4,nsl);
    pp->Iarray(pp->wet,nsl);
    pp->Iarray(pp->wet_n,nsl);
    pp->Iarray(pp->deep,nsl);
    for(int q=0; q<nsl; ++q)
    {
        pp->flagslice1[q] = pp->flagslice2[q] = pp->flagslice4[q] = 1;
        pp->wet[q] = pp->wet_n[q] = pp->deep[q] = 1;
    }

    // no boundary lists, no partition neighbours
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

    pp->vec2Dlength = nsl;
    pp->veclength = 0;
    pp->cellnum2D = c.nx*c.ny;
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

// fine bed: limited-linear prolongation of the parent bed (children average = parent)
void sflow_amr::ini_bed(lexer *p, sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm2D *pb = pfdm(c);
    fdm2D *b = c.b;

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        int ri=ii-1, rj=jj-1;
        int ic = c.pi0 + fdiv2(ri);
        int jc = c.pj0 + fdiv2(rj);
        double ox = (ri-2*fdiv2(ri)==0) ? -0.25 : 0.25;
        double oy = (rj-2*fdiv2(rj)==0) ? -0.25 : 0.25;

        double bc = pb->bed(ic,jc);
        double sx = mmod(pb->bed(ic+1,jc)-bc, bc-pb->bed(ic-1,jc));
        double sy = mmod(pb->bed(ic,jc+1)-bc, bc-pb->bed(ic,jc-1));

        b->bed(ii,jj) = bc + sx*ox + sy*oy;
        b->bed0(ii,jj) = b->bed(ii,jj);
        b->topobed(ii,jj) = b->bed(ii,jj);
        b->solidbed(ii,jj) = b->bed(ii,jj);
        b->depth(ii,jj) = p->wd - b->bed(ii,jj);
        b->ks(ii,jj) = p->B50;
    }

    // face depth as sflow_reconstruct::reconstruct_WL
    for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    b->dfx(ii,jj) = 0.5*(b->depth(ii+1,jj)+b->depth(ii,jj));

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
    b->dfy(ii,jj) = 0.5*(b->depth(ii,jj+1)+b->depth(ii,jj));
}

// --------------------------------------------------------------------- prolongation
// well balanced: the surface eta is interpolated with minmod slopes (switched off next
// to dry cells), the fine depth follows from the fine bed; momentum through u.
// Thin films (h_c < 3 A 244) copy the parent depth.
void sflow_amr::prolong(sflow_amr_patch &c, int ii, int jj, slice &WLp, slice &UHp, slice &VHp,
                        double &wl, double &uh, double &vh, int &w)
{
    fdm2D *pb = pfdm(c);
    const double wd = p0->A244;

    int ri=ii-1, rj=jj-1;
    int ic = c.pi0 + fdiv2(ri);
    int jc = c.pj0 + fdiv2(rj);
    double ox = (ri-2*fdiv2(ri)==0) ? -0.25 : 0.25;
    double oy = (rj-2*fdiv2(rj)==0) ? -0.25 : 0.25;

    auto get = [&](int a, int bb, double &e, double &u, double &v, int &wt)
    {
        double wlc = WLp(a,bb);
        e = wlc - pb->depth(a,bb);
        wt = pwet(c,a,bb);
        double wlvl = fabs(wlc)>wd ? wlc : 1.0e20;
        u = wt==1 ? UHp(a,bb)/wlvl : 0.0;
        v = wt==1 ? VHp(a,bb)/wlvl : 0.0;
    };

    double e0,u0,v0; int w0;
    get(ic,jc,e0,u0,v0,w0);

    if(w0!=1)
    {
        w=0; wl=wd; uh=vh=0.0;
        return;
    }

    double hc = WLp(ic,jc);
    if(hc<3.0*wd)
    {
        w=1; wl=hc; uh=UHp(ic,jc); vh=VHp(ic,jc);
        return;
    }

    double eE,uE,vE,eW,uW,vW,eN,uN,vN,eS,uS,vS;
    int wE,wW,wN,wS;
    get(ic+1,jc,eE,uE,vE,wE);
    get(ic-1,jc,eW,uW,vW,wW);
    get(ic,jc+1,eN,uN,vN,wN);
    get(ic,jc-1,eS,uS,vS,wS);

    double sxe=0.0,sye=0.0,sxu=0.0,syu=0.0,sxv=0.0,syv=0.0;
    if(wE==1 && wW==1 && wN==1 && wS==1)
    {
        sxe = mmod(eE-e0,e0-eW); sye = mmod(eN-e0,e0-eS);
        sxu = mmod(uE-u0,u0-uW); syu = mmod(uN-u0,u0-uS);
        sxv = mmod(vE-v0,v0-vW); syv = mmod(vN-v0,v0-vS);
    }

    double ef = e0 + sxe*ox + sye*oy;
    double hf = ef + c.b->depth(ii,jj);

    if(hf<=wd+eps)
    {
        w=0; wl=wd; uh=vh=0.0;
        return;
    }

    w=1;
    wl=hf;
    uh=hf*(u0 + sxu*ox + syu*oy);
    vh=hf*(v0 + sxv*ox + syv*oy);
}

// initial patch state from the parent, conservative for fully wet parents
void sflow_amr::ini_state(lexer *p, ghostcell *pgc, sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm2D *pb = pfdm(c);
    fdm2D *b = c.b;

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        double wl,uh,vh; int w;
        prolong(c,ii,jj,pb->WL,pb->UH,pb->VH,wl,uh,vh,w);
        b->WL(ii,jj)=wl; b->UH(ii,jj)=uh; b->VH(ii,jj)=vh;
        pp->wet[lij(pp,ii,jj)]=w;
    }

    // exact mass and momentum of each fully wet parent cell
    for(int ic=c.pi0; ic<=c.pi1; ++ic)
    for(int jc=c.pj0; jc<=c.pj1; ++jc)
    {
        int i0=1+2*(ic-c.pi0), j0=1+2*(jc-c.pj0);
        int nw = pp->wet[lij(pp,i0,j0)] + pp->wet[lij(pp,i0+1,j0)] + pp->wet[lij(pp,i0,j0+1)] + pp->wet[lij(pp,i0+1,j0+1)];
        if(nw<4)
        continue;

        double sw = b->WL(i0,j0)+b->WL(i0+1,j0)+b->WL(i0,j0+1)+b->WL(i0+1,j0+1);
        double su = b->UH(i0,j0)+b->UH(i0+1,j0)+b->UH(i0,j0+1)+b->UH(i0+1,j0+1);
        double sv = b->VH(i0,j0)+b->VH(i0+1,j0)+b->VH(i0,j0+1)+b->VH(i0+1,j0+1);
        double du = 4.0*pb->UH(ic,jc) - su;
        double dv = 4.0*pb->VH(ic,jc) - sv;
        double fw = 4.0*pb->WL(ic,jc)/sw;

        for(int a=0;a<2;++a)
        for(int d=0;d<2;++d)
        {
            b->UH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*du;
            b->VH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dv;
            b->WL(i0+a,j0+d) *= fw;
        }
    }

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        int w = pp->wet[lij(pp,ii,jj)];
        double wlvl = fabs(b->WL(ii,jj))>p->A244 ? b->WL(ii,jj) : 1.0e20;
        b->eta(ii,jj) = b->WL(ii,jj) - b->depth(ii,jj);
        b->eta_n(ii,jj) = b->eta(ii,jj);
        b->U(ii,jj) = w==1 ? b->UH(ii,jj)/wlvl : 0.0;
        b->V(ii,jj) = w==1 ? b->VH(ii,jj)/wlvl : 0.0;
        b->W(ii,jj) = 0.0;
        b->WH(ii,jj) = 0.0;
        b->hp(ii,jj) = b->WL(ii,jj);
        b->breaking(ii,jj) = 0;
        pp->wet_n[lij(pp,ii,jj)] = w;
    }

    // conserved variables and ghost cells as for level 0
    c.pmom->ini(pp,b,pgc);
}

// --------------------------------------------------------------------- stages
void sflow_amr::fill_ghosts(lexer *p, sflow_amr_patch &c, int s)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;
    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    slice *pWLi,*pUHi,*pVHi,*pWLo,*pUHo,*pVHo;

    c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);
    pmom(c)->stage_io(s,pfdm(c),pWLi,pUHi,pVHi,pWLo,pUHo,pVHo);

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        if(ii>=1 && ii<=c.nx && jj>=1 && jj<=c.ny)
        continue;

        double wl,uh,vh; int w;
        prolong(c,ii,jj,*pWLi,*pUHi,*pVHi,wl,uh,vh,w);

        (*WLi)(ii,jj)=wl; (*UHi)(ii,jj)=uh; (*VHi)(ii,jj)=vh;
        pp->wet[lij(pp,ii,jj)]=w;

        double wlvl = fabs(wl)>p->A244 ? wl : 1.0e20;
        b->eta(ii,jj) = wl - b->depth(ii,jj);
        b->U(ii,jj) = w==1 ? uh/wlvl : 0.0;
        b->V(ii,jj) = w==1 ? vh/wlvl : 0.0;
        b->W(ii,jj) = 0.0;
    }
}

void sflow_amr::restrict_patch(lexer *p, sflow_amr_patch &c, int s)
{
    fdm2D *b = c.b;
    fdm2D *pb = pfdm(c);
    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    slice *pWLi,*pUHi,*pVHi,*pWLo,*pUHo,*pVHo;

    c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);
    pmom(c)->stage_io(s,pb,pWLi,pUHi,pVHi,pWLo,pUHo,pVHo);

    for(int ic=c.pi0; ic<=c.pi1; ++ic)
    for(int jc=c.pj0; jc<=c.pj1; ++jc)
    {
        int i0=1+2*(ic-c.pi0), j0=1+2*(jc-c.pj0);

        double wl = 0.25*((*WLo)(i0,j0)+(*WLo)(i0+1,j0)+(*WLo)(i0,j0+1)+(*WLo)(i0+1,j0+1));
        double uh = 0.25*((*UHo)(i0,j0)+(*UHo)(i0+1,j0)+(*UHo)(i0,j0+1)+(*UHo)(i0+1,j0+1));
        double vh = 0.25*((*VHo)(i0,j0)+(*VHo)(i0+1,j0)+(*VHo)(i0,j0+1)+(*VHo)(i0+1,j0+1));

        (*pWLo)(ic,jc)=wl; (*pUHo)(ic,jc)=uh; (*pVHo)(ic,jc)=vh;

        int w = wl>p->A244+eps ? 1 : 0;
        set_pwet(c,ic,jc,w);

        double wlvl = fabs(wl)>p->A244 ? wl : 1.0e20;
        pb->eta(ic,jc) = wl - pb->depth(ic,jc);
        pb->U(ic,jc) = w==1 ? uh/wlvl : 0.0;
        pb->V(ic,jc) = w==1 ? vh/wlvl : 0.0;
        pb->W(ic,jc) = 0.0;
        pb->WH(ic,jc) = 0.0;
        pb->hp(ic,jc) = wl;
    }
}

void sflow_amr::step_begin(lexer *p, fdm2D *b, ghostcell *pgc)
{
    comms_off guard(pgc);

    for(auto &c : P)
    {
        c.pp->dt = p->dt;
        c.pp->dt_old = p->dt_old;
        c.pp->simtime = p->simtime;
        c.pp->count = p->count;
        c.pmom->inflow(c.pp,c.b,pgc,pflow_void);
    }
}

void sflow_amr::stage_begin(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(P.empty())
    return;

    comms_off guard(pgc);

    for(int n : order)
    fill_ghosts(p,P[n],s);

    // finest first: the parents' interface faces take the recorded fine fluxes
    for(int k=(int)order.size()-1; k>=0; --k)
    {
        sflow_amr_patch &c = P[order[k]];
        c.pmom->rk_stage(c.pp,c.b,pgc,s);
    }
}

void sflow_amr::stage_end(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(P.empty())
    return;

    comms_off guard(pgc);

    for(int k=(int)order.size()-1; k>=0; --k)
    restrict_patch(p,P[order[k]],s);
}

void sflow_amr::step_end(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(P.empty())
    return;

    comms_off guard(pgc);

    for(int n : order)
    P[n].pmom->rk_finish(P[n].pp,P[n].b,pgc);
}

// --------------------------------------------------------------------- flux matching
void sflow_amr::hll_hook(lexer *p, fdm2D *b, int ipol, int id)
{
    slice &Fx = (ipol==4) ? b->FEx : b->Fx;
    slice &Fy = (ipol==4) ? b->FEy : b->Fy;

    // a patch records its boundary faces
    if(id>0)
    {
        sflow_amr_patch &c = P[id-1];
        const int ie = c.nx, je = c.ny;   // last real cell: its high face is the boundary

        for(int r=0; r<c.ny; ++r)
        {
            c.rec[ipol][0][r] = Fx(0,r+1);
            c.rec[ipol][1][r] = Fx(ie,r+1);
        }
        for(int r=0; r<c.nx; ++r)
        {
            c.rec[ipol][2][r] = Fy(r+1,0);
            c.rec[ipol][3][r] = Fy(r+1,je);
        }

        if(ipol==4)
        {
            for(int r=0; r<c.ny; ++r)
            {
                c.rec[0][0][r] = b->dfx(0,r+1);
                c.rec[0][1][r] = b->dfx(ie,r+1);
            }
            for(int r=0; r<c.nx; ++r)
            {
                c.rec[0][2][r] = b->dfy(r+1,0);
                c.rec[0][3][r] = b->dfy(r+1,je);
            }
        }
    }

    // a parent takes the mean of the two fine faces on each coarse face of a child
    for(auto &c : P)
    {
        if(c.parent != id-1)
        continue;

        for(int jc=c.pj0; jc<=c.pj1; ++jc)
        {
            int r = 2*(jc-c.pj0);
            Fx(c.pi0-1,jc) = 0.5*(c.rec[ipol][0][r]+c.rec[ipol][0][r+1]);
            Fx(c.pi1,jc)   = 0.5*(c.rec[ipol][1][r]+c.rec[ipol][1][r+1]);

            if(ipol==4)
            {
                b->dfx(c.pi0-1,jc) = 0.5*(c.rec[0][0][r]+c.rec[0][0][r+1]);
                b->dfx(c.pi1,jc)   = 0.5*(c.rec[0][1][r]+c.rec[0][1][r+1]);
            }
        }
        for(int ic=c.pi0; ic<=c.pi1; ++ic)
        {
            int r = 2*(ic-c.pi0);
            Fy(ic,c.pj0-1) = 0.5*(c.rec[ipol][2][r]+c.rec[ipol][2][r+1]);
            Fy(ic,c.pj1)   = 0.5*(c.rec[ipol][3][r]+c.rec[ipol][3][r+1]);

            if(ipol==4)
            {
                b->dfy(ic,c.pj0-1) = 0.5*(c.rec[0][2][r]+c.rec[0][2][r+1]);
                b->dfy(ic,c.pj1)   = 0.5*(c.rec[0][3][r]+c.rec[0][3][r+1]);
            }
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
        for(auto &c : P)
        r = MIN(r, c.pp->DXM/p->DXM);

        r = pgc->globalmin(r);
        p->dt *= r;
        p->dt_old = p->dt;
        return;
    }

    // same CFL as sflow_etimestep, over the wet real patch cells
    const double g = fabs(p->W22);
    double cmin = 1.0e20;

    for(auto &c : P)
    {
        lexer *pp = c.pp;
        fdm2D *pb = c.b;
        for(int ii=1; ii<=c.nx; ++ii)
        for(int jj=1; jj<=c.ny; ++jj)
        if(pp->wet[lij(pp,ii,jj)]==1)
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

    if(p->mpirank==0 && (p->count%p->P12==0 || doprint))
    logout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<p->dt<<" \t "<<patches_total<<" \t "
          <<setprecision(15)<<m<<" \t "<<setprecision(6)<<(m-m0)/(fabs(m0)>1.0e-20?m0:1.0)<<endl;

    if(!doprint)
    return;

    write_vtr0(p,b);

    for(auto &c : P)
    write_vtr(p,c);

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
            out<<setprecision(10)<<f(ii,jj)+shift<<" ";
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

void sflow_amr::write_vtr(lexer *p, sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm2D *pb = c.b;
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i-%04i-%04i.vtr",printcount_amr,p->mpirank+1,c.id+1);

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
        for(int jj=1; jj<=c.ny; ++jj)
        {
            for(int ii=1; ii<=c.nx; ++ii)
            out<<setprecision(10)<<f(ii,jj)+shift<<" ";
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

    out<<"<DataArray type=\"Int32\" Name=\"wetdry\" format=\"ascii\">\n";
    for(int jj=1; jj<=c.ny; ++jj)
    {
        for(int ii=1; ii<=c.nx; ++ii)
        out<<pp->wet[lij(pp,ii,jj)]<<" ";
        out<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"</CellData>\n<Coordinates>\n";

    const int m = marge;
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=1; ii<=c.nx+1; ++ii) out<<setprecision(12)<<pp->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=1; jj<=c.ny+1; ++jj) out<<setprecision(12)<<pp->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}
