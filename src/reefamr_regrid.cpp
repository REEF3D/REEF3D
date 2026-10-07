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

#include"reefamr.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<cmath>
#include<iomanip>
#include<algorithm>
#include<iostream>

namespace
{
// floor(a/2^l) for negative a as well
inline int fsh(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }
}

//  Regridding: the patches are rebuilt from the refinement flags of the module (tag), the
//  refinement boxes and the zone of the moving bodies, with a buffer of nbuf cells.  Flags
//  mark tiles of the global index space of each level; the tile maps are global, so the
//  refined region does not depend on the domain decomposition.  Marked tiles are merged into
//  rectangles and cut at the partition edges.  A patch with the same box as before is kept
//  with its state.

void reefamr::regrid(lexer *p, ghostcell *pgc, bool initial)
{
    const int T = tile;
    const int old_total = patches_total;

    // current hierarchy: final state with filled ghost cells
    regrid_prepare(pgc);

    // ---- tile maps, finest first: flags of the next coarser level, the refinement boxes
    //      and the footprint of the finer level grown by the nesting width
    vector<vector<unsigned char>> M(maxlev+1);

    if(par.zones)
    zone_setup(p);

    // cells that no patch may cover now (dry or shallow cells of the module), grown by dryband
    const bool dyn = par.dryband>0;
    if(dyn)
    {
        const int B = par.dryband;
        forbidd.assign((size_t)GNX*GNY,0);
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(cell_unfit(ii,jj))
        {
            int I=ii+O0i, J=jj+O0j;
            for(int a=MAX(I-B,0); a<=MIN(I+B,GNX-1); ++a)
            for(int d=MAX(J-B,0); d<=MIN(J+B,GNY-1); ++d)
            forbidd[(size_t)a*GNY+d]=1;
        }
        global_or(forbidd);
    }

    auto tile_forbidden = [&](int l, int ti, int tj)
    {
        int a0 = (ti*T)>>l, a1 = (MIN((ti+1)*T,GNX<<l)-1)>>l;
        int d0 = (tj*T)>>l, d1 = (MIN((tj+1)*T,GNY<<l)-1)>>l;
        for(int a=a0; a<=a1; ++a)
        for(int d=d0; d<=d1; ++d)
        if(forbid0[(size_t)a*GNY+d] || (dyn && forbidd[(size_t)a*GNY+d]))
        return true;
        return false;
    };

    for(int l=maxlev; l>=1; --l)
    {
        M[l].assign((size_t)gtnx[l]*gtny[l],0);

        tag(l-1,M[l]);

        for(size_t q=0; q+3<par.rbox.size(); q+=4)
        SLICELOOP4
        if(p->XP[IP]>=par.rbox[q] && p->XP[IP]<=par.rbox[q+1] && p->YP[JP]>=par.rbox[q+2] && p->YP[JP]<=par.rbox[q+3])
        {
            int I=i+O0i, J=j+O0j;
            for(int ti=(I<<l)/T; ti<=(((I+1)<<l)-1)/T; ++ti)
            for(int tj=(J<<l)/T; tj<=(((J+1)<<l)-1)/T; ++tj)
            M[l][(size_t)ti*gtny[l]+tj]=1;
        }

        if(par.zones)
        SLICELOOP4
        if(zone_test(p->XP[IP],p->YP[JP]))
        {
            int I=i+O0i, J=j+O0j;
            for(int ti=(I<<l)/T; ti<=(((I+1)<<l)-1)/T; ++ti)
            for(int tj=(J<<l)/T; tj<=(((J+1)<<l)-1)/T; ++tj)
            M[l][(size_t)ti*gtny[l]+tj]=1;
        }

        global_or(M[l]);

        // hysteresis: a refined tile is kept until it has not been flagged for keep
        // regrids, so that a flag which drops below the threshold for a moment does not remove
        // and rebuild the patch
        if(keep>0)
        {
            vector<unsigned char> &A = tage[l];
            if(A.size()!=M[l].size())
            A.assign(M[l].size(),255);

            const bool old = gtile[l].size()==M[l].size();
            for(size_t t=0; t<M[l].size(); ++t)
            {
                if(M[l][t])
                A[t] = 0;
                else if(A[t]<255)
                ++A[t];

                if(!initial && old && !M[l][t] && gtile[l][t] && A[t]<=keep)
                M[l][t] = 1;
            }
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

        for(int ti=0; ti<gtnx[l]; ++ti)
        for(int tj=0; tj<gtny[l]; ++tj)
        if(M[l][(size_t)ti*gtny[l]+tj] && tile_forbidden(l,ti,tj))
        M[l][(size_t)ti*gtny[l]+tj]=0;

        // lazy layout: the refined tiles of the current layout stay as long as the union with the
        // flagged tiles is at most par.lazy times the flagged tiles.  A patch is kept only with
        // exactly its old box, so every tile that comes or goes rebuilds the patches of the rank
        // (lexer, fdm, kernels, multigrid, body grids); with the union a zone that oscillates
        // with the body settles after one period.  The map is global, all ranks decide alike
        if(par.lazy>0.0 && !initial && gtile[l].size()==M[l].size())
        {
            long nflag=0, nunion=0;
            for(size_t t=0; t<M[l].size(); ++t)
            {
                nflag += M[l][t] ? 1 : 0;
                nunion += (M[l][t] || gtile[l][t]) ? 1 : 0;
            }

            // (a tile of the layout that is forbidden now, e.g. fallen dry, is not kept)
            if(nflag>0 && double(nunion)<=par.lazy*double(nflag))
            for(size_t t=0; t<M[l].size(); ++t)
            if(gtile[l][t] && !(dyn && tile_forbidden(l,(int)(t/gtny[l]),(int)(t%gtny[l]))))
            M[l][t] = 1;
        }
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

    // ---- lazy layout: the tile maps did not change, the hierarchy stays as it is
    if(par.lazy>0.0 && !initial)
    {
        bool same = true;
        for(int l=1; l<=maxlev && same; ++l)
        if(gtile[l]!=M[l])
        same = false;

        if(same)
        {
            ++regrids_skipped;
            return;
        }
    }

    // ---- patches: marked tiles merged into rectangles, cut at the rank box
    vector<reefamr_patch*> oldP = P;
    vector<char> kept(oldP.size(),0);
    vector<reefamr_patch*> newP;
    vector<vector<int>> newlev(maxlev+1);

    if(par.place>0)
    place_patches(p,pgc,M,oldP,kept,newP,newlev);
    else
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

            reefamr_patch *c=nullptr;
            for(size_t k=0; k<oldP.size(); ++k)
            if(!kept[k] && oldP[k]->lev==l && oldP[k]->I0==I0 && oldP[k]->I1==I1 && oldP[k]->J0==J0 && oldP[k]->J1==J1)
            {
                kept[k]=1;
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

    vector<reefamr_patch*> gone;
    for(size_t k=0; k<oldP.size(); ++k)
    if(!kept[k])
    gone.push_back(oldP[k]);

    P = newP;
    lev = newlev;
    regrid_ids();

    build_tiles();
    build_gtable();

    // ---- time-independent data of the new patches (bed)
    regrid_static(pgc);

    // ---- fill, restriction and flux matching plans of the new hierarchy
    build_plans(pgc);

    // ---- state of the new patches, coarse to fine
    regrid_state(pgc,oldP);

    for(auto c : gone)
    free_patch(c);

    for(auto c : P)
    c->fresh = false;

    long cells=0;
    for(auto c : P)
    cells += (long)c->nx*c->ny;

    patches_total = pgc->globalisum((int)P.size());
    nlevg.assign(maxlev+1,0);
    for(int l=1; l<=maxlev; ++l)
    nlevg[l] = pgc->globalisum((int)lev[l].size());
    cells_total = (long)pgc->globalsum(double(cells));
    cells_local = cells;

    // ---- level 0 consistent with the patches
    regrid_finish(pgc,old_total);

    ++regrids;
}

// --------------------------------------------------------------------- plans
void reefamr::build_plans(ghostcell *pgc)
{
    const int me = p0->mpirank;
    const int np = p0->mpi_size;

    match.assign(P.size()+1,vector<reefamr_match>());
    rmatch.assign(P.size()+1,vector<reefamr_match>());

    // local source of a cell (l,I,J) held by this rank: a patch of level l or the coarser level
    int violations=0;
    auto source = [&](int l, int I, int J, reefamr_fill &f)
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
        int oi,oj;
        goff(g,oi,oj);
        f.kind = 1;
        f.g = g;
        f.si = Ic-oi;
        f.sj = Jc-oj;
        f.ox = (I-2*Ic==0) ? -1 : 1;
        f.oy = (J-2*Jc==0) ? -1 : 1;
        return true;
    };

    for(int l=1; l<=maxlev; ++l)
    {
        // ---- cells around the patches
        vector<vector<int>> req(np), srv;
        vector<vector<int>> dst(np);

        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            lexer *pp = c->pp;
            c->fill.clear();

            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                if(ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny)
                continue;

                if(pp->flagslice4[(ii-pp->imin)*pp->jmax + (jj-pp->jmin)]<0)
                continue;

                int I = ii-EXT+c->I0, J = jj-EXT+c->J0;
                int o = holder(l,I,J);
                if(o<0)
                continue;

                reefamr_fill f;
                f.di=ii; f.dj=jj; f.g=-1; f.si=f.sj=0; f.ox=f.oy=0; f.slot=-1;
                f.aux = fill_aux(c,ii,jj);

                if(o==me)
                {
                    if(source(l,I,J,f))
                    c->fill.push_back(f);
                    continue;
                }

                f.kind = 2;
                c->fill.push_back(f);
                req[o].push_back(I); req[o].push_back(J);
                dst[o].push_back(id); dst[o].push_back((int)c->fill.size()-1);
            }
        }

        xsetup(req,srv,2);

        reefamr_xplan &X = gplan[l];
        X = reefamr_xplan();
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
                reefamr_fill f;
                f.di=f.dj=0; f.g=-1; f.si=f.sj=0; f.ox=f.oy=0; f.slot=-1;
                f.aux = serve_aux(l,I,J);
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
                reefamr_fill *f = &P[dst[r][k]]->fill[dst[r][k+1]];
                f->si = kp;
                f->slot = (int)(k/2);
                grecv[l].push_back(f);
            }
        }

        // ---- restriction targets: the parent cell of every 2x2 block, on this rank or (-3) on
        //      the rank that holds it, through the block plans
        vector<vector<int>> breq(np), bsrq, bdst(np);
        bnloc[l] = 0;

        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            int nb = (c->nx/2)*(c->ny/2);
            c->rgrid.assign(nb,-2); c->ric.assign(nb,0); c->rjc.assign(nb,0);

            for(int bi=0; bi<c->nx/2; ++bi)
            for(int bj=0; bj<c->ny/2; ++bj)
            {
                int k = bi*(c->ny/2)+bj;
                int Ic = (c->I0>>1)+bi, Jc = (c->J0>>1)+bj;
                if(flag0(fsh(Ic,l-1),fsh(Jc,l-1))<0)
                continue;
                int o = holder(l-1,Ic,Jc);
                if(o>=0 && o!=me)
                {
                    c->rgrid[k]=-3;
                    breq[o].push_back(Ic); breq[o].push_back(Jc);
                    bdst[o].push_back(id); bdst[o].push_back(k);
                    continue;
                }
                int g = patch_at(l-1,Ic,Jc);
                if(g<-1)
                {
                    ++violations;
                    continue;
                }
                int oi,oj;
                goff(g,oi,oj);
                c->rgrid[k]=g;
                c->ric[k]=Ic-oi;
                c->rjc[k]=Jc-oj;
                ++bnloc[l];
            }
        }

        xsetup(breq,bsrq,2);

        reefamr_xplan &BU = bup[l];
        reefamr_xplan &BD = bdn[l];
        BU = reefamr_xplan();
        BD = reefamr_xplan();
        bloc[l].clear();
        bsrv[l].clear();

        // my blocks with a parent on rank r: up to r, down from r
        for(int r=0; r<np; ++r)
        if(!breq[r].empty())
        {
            vector<int> items;
            for(size_t k=0; k<bdst[r].size(); k+=2)
            {
                items.push_back((int)bloc[l].size()/2);
                bloc[l].push_back(bdst[r][k]);
                bloc[l].push_back(bdst[r][k+1]);
            }
            BU.speer.push_back(r);
            BU.sitem.push_back(items);
            BU.sbuf.push_back(vector<double>());
            BD.rpeer.push_back(r);
            BD.rcount.push_back((int)breq[r].size()/2);
        }

        // the parents I hold for blocks of rank r
        for(int r=0; r<np; ++r)
        if(!bsrq[r].empty())
        {
            vector<int> items;
            for(size_t k=0; k<bsrq[r].size(); k+=2)
            {
                int Ic=bsrq[r][k], Jc=bsrq[r][k+1];
                reefamr_block b{-2,0,0};
                int g = patch_at(l-1,Ic,Jc);
                if(g>=-1)
                {
                    int oi,oj;
                    goff(g,oi,oj);
                    b = reefamr_block{g,Ic-oi,Jc-oj};
                }
                else
                ++violations;
                items.push_back((int)bsrv[l].size());
                bsrv[l].push_back(b);
            }
            BU.rpeer.push_back(r);
            BU.rcount.push_back((int)bsrq[r].size()/2);
            BD.speer.push_back(r);
            BD.sitem.push_back(items);
            BD.sbuf.push_back(vector<double>());
        }

        // ---- flux matching: coarse faces next to the patch boundary
        vector<vector<int>> freq(np), fsrv;
        fsend[l].clear();
        vector<vector<int>> fitems(np);

        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];

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

                    int o = holder(l-1,Ic,Jc);
                    if(o<0)
                    continue;

                    if(o==me)
                    {
                        if(covered(l,2*Ic,2*Jc))
                        continue;

                        int g = patch_at(l-1,Ic,Jc);
                        if(g<-1)
                        {
                            ++violations;
                            continue;
                        }
                        int oi,oj;
                        goff(g,oi,oj);
                        reefamr_match mt;
                        mt.dir = side<2 ? 0 : 1;
                        mt.fi = Ic-oi - (side==1 ? 1 : 0);
                        mt.fj = Jc-oj - (side==3 ? 1 : 0);
                        mt.child = id;
                        mt.side = side;
                        mt.r = 2*k;
                        for(int q=0;q<reefamr_match::NVAL;++q) mt.val[q]=0.0;
                        match[g+1].push_back(mt);
                        continue;
                    }

                    freq[o].push_back(side); freq[o].push_back(Ic); freq[o].push_back(Jc);
                    fitems[o].push_back(id); fitems[o].push_back(side); fitems[o].push_back(2*k);
                }
            }
        }

        xsetup(freq,fsrv,3);

        reefamr_xplan &F = fplan[l];
        F = reefamr_xplan();
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

                if(!covered(l,2*Ic,2*Jc))
                {
                    int g = patch_at(l-1,Ic,Jc);
                    if(g>=-1)
                    {
                        int oi,oj;
                        goff(g,oi,oj);
                        reefamr_match mt;
                        mt.dir = side<2 ? 0 : 1;
                        mt.fi = Ic-oi - (side==1 ? 1 : 0);
                        mt.fj = Jc-oj - (side==3 ? 1 : 0);
                        mt.child = -1;
                        mt.side = side;
                        mt.r = 0;
                        for(int q=0;q<reefamr_match::NVAL;++q) mt.val[q]=0.0;
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
    cout<<par.name<<": "<<violations<<" cells without a coarser grid (nesting)"<<endl;
}

// --------------------------------------------------------------------- placement (par.place 1)
//  The marked tiles of every level are merged into rectangles over the whole domain (not cut at
//  the rank boxes), the same on every rank.  Rectangles whose work is more than a quarter of the
//  mean work of a rank are split along their longer side at tile edges (finer pieces balance
//  better on few ranks, at the price of more patch edges).  The work of a rank is its level-0
//  cells plus the cells of its patches (EXT cells and layers included).  The pieces are placed
//  largest first: on the rank of the level-0 cell at their centre if it stays within 5 % of the
//  mean work, else on the rank with the least work.  The previous placement is kept (kept boxes on
//  their rank, new pieces placed as above) unless the new one lowers the predicted maximum work by
//  more than par.rebalance.  Every rank computes the same placement and makes its own patches.
void reefamr::place_patches(lexer *p, ghostcell *pgc, vector<vector<unsigned char>> &M, vector<reefamr_patch*> &oldP,
                            vector<char> &kept, vector<reefamr_patch*> &newP, vector<vector<int>> &newlev)
{
    const int T = tile;
    const int me = p->mpirank;
    const int np = p->mpi_size;

    struct piece { int lev, I0, I1, J0, J1; double w; int rank; };
    vector<piece> pcs;

    auto work = [&](int l, int I0, int I1, int J0, int J1)
    {
        return double(I1-I0+1+2*EXT)*double(J1-J0+1+2*EXT)*double(p->knoz*vfac[l]);
    };

    // rectangles of marked tiles, all levels
    for(int l=1; l<=maxlev; ++l)
    {
        gtile[l] = M[l];

        struct R { int ta,tb,tj0,tj1; };
        vector<R> open, done;

        for(int tj=0; tj<gtny[l]; ++tj)
        {
            vector<R> next;
            int ti=0;
            while(ti<gtnx[l])
            {
                if(!M[l][(size_t)ti*gtny[l]+tj]) { ++ti; continue; }
                int ta=ti;
                while(ti<gtnx[l] && M[l][(size_t)ti*gtny[l]+tj]) ++ti;
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
        pcs.push_back({l,o.ta,o.tb,o.tj0,o.tj1,0.0,-1});      // tile ranges for now
    }

    // work per rank on level 0 and the mean work per rank
    vector<double> base(np);
    double wtot = 0.0;
    for(int r=0; r<np; ++r)
    {
        base[r] = double(rbx1[r]-rbx0[r]+1)*double(rby1[r]-rby0[r]+1)*double(p->knoz);
        wtot += base[r];
    }
    auto box = [&](const piece &c, int &I0, int &I1, int &J0, int &J1)
    {
        I0 = c.I0*T; I1 = MIN((c.I1+1)*T-1,(GNX<<c.lev)-1);
        J0 = c.J0*T; J1 = MIN((c.J1+1)*T-1,(GNY<<c.lev)-1);
    };
    for(auto &c : pcs)
    {
        int I0,I1,J0,J1;
        box(c,I0,I1,J0,J1);
        wtot += work(c.lev,I0,I1,J0,J1);
    }
    const double wmean = wtot/np;
    const double wmax = (par.place==2 ? 0.125 : 0.25)*wmean;

    // split the large rectangles at tile edges, then cell boxes, only pieces with fluid
    vector<piece> parts;
    for(auto &c : pcs)
    {
        int I0,I1,J0,J1;
        box(c,I0,I1,J0,J1);
        const double w = work(c.lev,I0,I1,J0,J1);
        const int nti = c.I1-c.I0+1, ntj = c.J1-c.J0+1;
        int n = MAX(1,(int)ceil(w/wmax));

        // rectangles at a body: whole (par.place_whole_zones)
        if(par.place_whole_zones)
        for(const reefamr_zone &z : zones)
        {
            const double x0 = gnode(0,c.lev,I0), x1 = gnode(0,c.lev,I1+1);
            const double y0 = gnode(1,c.lev,J0), y1 = gnode(1,c.lev,J1+1);
            if(x1>=z.bx0-par.zr && x0<=z.bx1+par.zr && (p->j_dir==0 || (y1>=z.by0-par.zr && y0<=z.by1+par.zr)))
            n = 1;
        }
        const bool alongx = (nti>=ntj);
        n = MIN(n,alongx ? nti : ntj);

        for(int s=0; s<n; ++s)
        {
            piece q = c;
            if(alongx)
            {
                q.I0 = c.I0 + (s*nti)/n;
                q.I1 = c.I0 + ((s+1)*nti)/n - 1;
            }
            else
            {
                q.J0 = c.J0 + (s*ntj)/n;
                q.J1 = c.J0 + ((s+1)*ntj)/n - 1;
            }
            int a0,a1,b0,b1;
            box(q,a0,a1,b0,b1);

            int fluid=0;
            for(int a=(a0>>q.lev); a<=(a1>>q.lev) && fluid==0; ++a)
            for(int d=(b0>>q.lev); d<=(b1>>q.lev) && fluid==0; ++d)
            if(flag0(a,d)>0)
            fluid=1;
            if(fluid==0)
            continue;

            q.I0=a0; q.I1=a1; q.J0=b0; q.J1=b1;
            q.w = work(q.lev,a0,a1,b0,b1);
            parts.push_back(q);
        }
    }

    // the mean work per rank of the pieces (their EXT cells add to the work of the rectangles)
    double wpl = 0.0;
    for(int r=0; r<np; ++r)
    wpl += base[r];
    for(auto &q : parts)
    wpl += q.w;
    const double wmeanp = wpl/np;

    // largest first (ties: level, box), deterministic
    vector<int> ord(parts.size());
    for(size_t k=0; k<ord.size(); ++k)
    ord[k]=(int)k;
    std::sort(ord.begin(),ord.end(),[&](int a, int b)
    {
        const piece &A = parts[a], &B = parts[b];
        if(A.w!=B.w) return A.w>B.w;
        if(A.lev!=B.lev) return A.lev<B.lev;
        if(A.I0!=B.I0) return A.I0<B.I0;
        return A.J0<B.J0;
    });

    auto least = [&](const vector<double> &ld)
    {
        int r=0;
        for(int q=1; q<np; ++q)
        if(ld[q]<ld[r])
        r=q;
        return r;
    };

    auto greedy = [&](const piece &c, vector<double> &ld)
    {
        const int o = owner(((c.I0+c.I1)/2)>>c.lev,((c.J0+c.J1)/2)>>c.lev);
        if(par.place==1 && o>=0 && ld[o]+c.w<=1.05*wmeanp)
        return o;
        return least(ld);
    };

    // the new placement
    vector<int> rnew(parts.size());
    vector<double> lnew = base;
    for(int k : ord)
    {
        rnew[k] = greedy(parts[k],lnew);
        lnew[rnew[k]] += parts[k].w;
    }

    // the previous placement: boxes of the last layout keep their rank
    vector<int> rold(parts.size(),-1);
    vector<double> lold = base;
    bool anyold = false;
    for(int k : ord)
    {
        const piece &c = parts[k];
        for(const reefamr_gpatch &G : GP)
        if(G.lev==c.lev && G.I0==c.I0 && G.I1==c.I1 && G.J0==c.J0 && G.J1==c.J1)
        {
            rold[k] = G.rank;
            anyold = true;
            break;
        }
        if(rold[k]>=0)
        lold[rold[k]] += c.w;
    }
    for(int k : ord)
    if(rold[k]<0)
    {
        rold[k] = greedy(parts[k],lold);
        lold[rold[k]] += parts[k].w;
    }

    double mnew=0.0, mold=0.0;
    for(int r=0; r<np; ++r)
    {
        mnew = MAX(mnew,lnew[r]);
        mold = MAX(mold,lold[r]);
    }
    const bool fresh = !anyold || (mold-mnew > par.rebalance*mold) || par.place==2;
    const vector<int> &rk = fresh ? rnew : rold;

    // my patches, level by level in the order of the pieces (kept: same box, on this rank before)
    for(int l=1; l<=maxlev; ++l)
    for(size_t k=0; k<parts.size(); ++k)
    {
        const piece &c = parts[k];
        if(c.lev!=l || rk[k]!=me)
        continue;

        reefamr_patch *q=nullptr;
        for(size_t m=0; m<oldP.size(); ++m)
        if(!kept[m] && oldP[m]->lev==l && oldP[m]->I0==c.I0 && oldP[m]->I1==c.I1 && oldP[m]->J0==c.J0 && oldP[m]->J1==c.J1)
        {
            kept[m]=1;
            q = oldP[m];
            q->fresh = false;
            break;
        }

        if(q==nullptr)
        q = make_patch(p,pgc,l,c.I0,c.I1,c.J0,c.J1);

        newlev[l].push_back((int)newP.size());
        newP.push_back(q);
    }

    if(me==0 && regrids>0 && fresh && anyold)
    cout<<par.name<<": patches placed anew, predicted load max/mean "<<setprecision(3)<<mold/wmeanp<<" -> "<<mnew/wmeanp<<endl;
}
