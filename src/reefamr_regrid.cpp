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

            if(nflag>0 && double(nunion)<=par.lazy*double(nflag))
            for(size_t t=0; t<M[l].size(); ++t)
            if(gtile[l][t])
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

        // ---- restriction targets
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
            }
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
                        if(patch_at(l,2*Ic,2*Jc)>=0)
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

                if(patch_at(l,2*Ic,2*Jc)<0)
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
