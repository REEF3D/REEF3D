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

//  SFLOW mesh refinement in checkpoints (P 48, I 45): the arrays of sflow_ckpt_fields and of the
//  momentum, free-surface, pressure and diffusion objects of every patch (reefamr_ckpt.cpp), the
//  hierarchy of the checkpoint restored by a regrid to its tile maps.
//  Sediment (S 10, sflow_amr_sed.cpp): the slices of the sediment objects of a patch (and of a block
//  of a multi-block grid) from the registry of its lexer, as the CFD patches (cfd_amr_ckpt.cpp), and
//  the scalars of its sediment_f (sediment_f::ckpt_state) per patch: from the patch of the checkpoint
//  with the same level and box, else from the one of the level that holds the centre of the patch.
//  The run scalars of the patch lexers (time, step, the sediment time and step, ...: the patches take
//  them over from level 0 during the step, some only after they are used) in the same way, and the output
//  counters of the refinement (the interval P 182 / P 181 of the log and the grids).

#include"sflow_amr.h"
#include"sflow_momentum_SI.h"
#include"reef3d_mpmd.h"
#include"reefckpt.h"
#include"sflow_ckpt_fields.h"
#include"sflow_momentum_RK3.h"
#include"sflow_fsf.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"
#include"sediment_f.h"
#include"ckpt_registry.h"
#include"lexer.h"
#include<set>
#include<cmath>
#include<mpi.h>

using namespace std;

namespace
{
// the slices and cell fields of the sediment objects of a patch: the registry of its lexer, without the
// arrays of the explicit list (kinds S, C; the solver arrays local to the rank are not kept)
template<class F>
void sed_registry(lexer *pp, const std::set<const void*> &listed, F &&f)
{
    ckpt_registry::visit(pp,[&](const ckpt_registry::entry &e)
    {
        if(e.local || (e.kind!='S' && e.kind!='C') || listed.count(e.V))
        return;
        f(e);
    });
}

// the run scalars of a patch lexer
template<class F>
void lexer_scalars(lexer *q, F &&f)
{
    f(q->simtime); f(q->dt); f(q->dt_old); f(q->sedtime); f(q->dtsed); f(q->sedsimtime); f(q->sedwavetime);
    f(q->wavetime); f(q->bedmin); f(q->bedmax); f(q->slidecells);
    double c = q->count;
    f(c);
    q->count = int(c);
}

void patch_fields(sflow_amr_patch *c, const reefamr::ckpt_visit &v)
{
    std::set<const void*> listed;
    auto fd = [&](const char *n, char k, double *V){ listed.insert(V); v(n,k,V,nullptr); };
    auto fi = [&](const char *n, char k, int *V){ listed.insert(V); v(n,k,nullptr,V); };
    sflow_ckpt_fields(c->pp,c->b,fd,fi);
    // shoreline: the restricted surface offsets of the covered cells (used by the stages of the next
    // step) and the parent surfaces of the ghost cells
    if(c->sfc.size()==(size_t)c->pp->imax*c->pp->jmax)
    fd("amr_sfc",'S',c->sfc.data());
    if(c->gsf.size()==(size_t)c->pp->imax*c->pp->jmax)
    fd("amr_gsf",'S',c->gsf.data());
    if(c->pmom!=nullptr)
    c->pmom->ckpt_fields(fd,fi);
    if(c->psim!=nullptr)
    c->psim->ckpt_fields(fd,fi);
    if(c->pfsf!=nullptr)
    c->pfsf->ckpt_fields(fd,fi);
    if(c->ppress!=nullptr)
    c->ppress->ckpt_fields(fd,fi);
    if(c->pdiff!=nullptr)
    c->pdiff->ckpt_fields(fd,fi);

    // S 10: the sediment cell flags and the arrays of the sediment objects of the patch
    if(c->sed!=nullptr)
    {
        if(c->pp->DFBED!=nullptr)
        fi("amr_DFBED",'S',c->pp->DFBED);
        sed_registry(c->pp,listed,[&](const ckpt_registry::entry &e)
        {
            if(e.isint) v("amr_sedreg",e.kind,nullptr,(int*)e.V);
            else        v("amr_sedreg",e.kind,(double*)e.V,nullptr);
        });
    }
}

// the code of the arrays of a patch (its count and kinds), to check that the patches of the run take the
// arrays of the checkpoint in the same order
double fields_code(sflow_amr_patch *c)
{
    double code = 0.0;
    int n = 0;
    patch_fields(c,[&](const char*, char k, double *V, int *I)
    {
        ++n;
        code = code*3.0 + (k=='S' ? 1.0 : 2.0) + (V==nullptr && I!=nullptr ? 0.5 : 0.0);
        code = std::fmod(code,1.0e12);
    });
    return code + 1.0e12*n;
}
}

void sflow_amr::ckpt_write(reefckpt &ck)
{
    ckpt_write_layout(ck);
    ckpt_write_patches(ck,"sflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(SP(q),v); });

    // the output counters (the same on all ranks)
    ck.blob("sflow_amr_print",vector<double>{printtime_amr,double(printcount_amr)});

    // the scalars of every patch of this rank: level, box, the code of its arrays, the run scalars of its
    // lexer, the length and the values of sediment_f::ckpt_state (S 10)
    vector<double> b;
    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        if(q->lev<1 && q->blk==0)
        continue;
        vector<double> lx;
        lexer_scalars(q->pp,[&](double &x){ lx.push_back(x); });
        vector<double> st;
        if(c->sed!=nullptr)
        c->sed->ckpt_state(c->pp,st);
        const double hd[6] = {double(q->lev),double(q->I0),double(q->I1),double(q->J0),double(q->J1),fields_code(c)};
        b.insert(b.end(),hd,hd+6);
        b.push_back(double(lx.size()));
        b.insert(b.end(),lx.begin(),lx.end());
        b.push_back(double(st.size()));
        b.insert(b.end(),st.begin(),st.end());
        b.push_back(double(c->fill.size()));
        for(const auto &f : c->fill)
        b.push_back(f.aux);
    }
    ck.rank_blob("sflow_amr_patches",b);

    // the still water depths of the filled cells this rank serves and of the ring of the native grid (multi-block):
    // with sediment they follow the bed (sed_fill_bed), the regrid of the restart takes them from the bed of the
    // coarser level
    vector<double> a;
    a.push_back(double(gserve.size()));
    for(const auto &G : gserve)
    {
        a.push_back(double(G.size()));
        for(const auto &f : G)
        a.push_back(f.aux);
    }
    a.push_back(double(nfill.size()));
    for(const auto &f : nfill)
    a.push_back(f.aux);
    ck.rank_blob("sflow_amr_aux",a);
}

bool sflow_amr::ckpt_read(reefckpt &ck, string &msg)
{
    // the hierarchy of the checkpoint (a regrid to its tile maps, the state of level 0 read before),
    // then the arrays of its patches
    if(!ckpt_read_layout(ck,msg))
    return false;
    regrid(p0,pgc0,false);

    if(!ckpt_read_patches(ck,"sflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(SP(q),v); },msg))
    return false;

    // semi-implicit free surface: the face discharges of the checkpoint
    for(auto q : P)
    if(SP(q)->psim!=nullptr)
    SP(q)->psim->ckpt_restored();

    if(ck.has("sflow_amr_print"))
    {
        const vector<double> pr = ck.get_blob("sflow_amr_print");
        if(pr.size()==2)
        {
            printtime_amr = pr[0];
            printcount_amr = int(pr[1]);
        }
    }

    // the scalars of the patches (lexer, sediment), the still water depths of their filled cells, and the
    // check of the arrays of the same patches
    struct sp { int lev, I0, I1, J0, J1; double code; vector<double> lx, st, fa; };
    vector<sp> cs;
    bool torn = false;
    for(const auto &bl : ck.get_rank_blobs("sflow_amr_patches"))
    {
        size_t k=0;
        while(k+6<=bl.size() && !torn)
        {
            sp e;
            e.lev = int(bl[k]); e.I0 = int(bl[k+1]); e.I1 = int(bl[k+2]); e.J0 = int(bl[k+3]); e.J1 = int(bl[k+4]);
            e.code = bl[k+5];
            k += 6;
            for(vector<double> *a : {&e.lx,&e.st,&e.fa})
            {
                const size_t n = k<bl.size() ? size_t(bl[k]) : 0;
                ++k;
                if(k+n>bl.size())
                {
                    torn = true;
                    break;
                }
                a->assign(bl.begin()+k,bl.begin()+k+n);
                k += n;
            }
            cs.push_back(e);
        }
    }
    if(torn)
    {
        msg = "the patches of the mesh refinement in the checkpoint are incomplete";
        return false;
    }

    int bad = 0;
    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        if(q->lev<1 && q->blk==0)
        continue;

        // the same patch, else the patch of the level that holds the centre of this one
        const sp *src = nullptr;
        bool same = false;
        for(const auto &e : cs)
        if(e.lev==q->lev && e.I0==q->I0 && e.I1==q->I1 && e.J0==q->J0 && e.J1==q->J1)
        {
            src = &e;
            same = true;
            if(e.code!=fields_code(c))
            bad = 1;
            break;
        }
        const int Ic = (q->I0+q->I1)/2, Jc = (q->J0+q->J1)/2;
        if(src==nullptr)
        for(const auto &e : cs)
        if(e.lev==q->lev && Ic>=e.I0 && Ic<=e.I1 && Jc>=e.J0 && Jc<=e.J1)
        {
            src = &e;
            break;
        }

        // the filled cells: their still water depth as in the checkpoint, else from the depth of the patch
        if(same && src->fa.size()==c->fill.size())
        for(size_t m=0; m<c->fill.size(); ++m)
        c->fill[m].aux = src->fa[m];
        else
        for(auto &f : c->fill)
        f.aux = c->b->depth(f.di,f.dj);

        if(src==nullptr)
        continue;

        size_t m = 0;
        lexer_scalars(q->pp,[&](double &x){ if(m<src->lx.size()) x = src->lx[m]; ++m; });
        if(c->sed!=nullptr && !c->sed->ckpt_restore(c->pp,pgc0,src->st))
        bad = 1;
    }
    int badall = 0;
    MPI_Allreduce(&bad,&badall,1,MPI_INT,MPI_MAX,reef3d_comm());
    if(badall)
    {
        msg = "the patches of the mesh refinement in the checkpoint do not fit this run (other keys?)";
        return false;
    }

    // the served cells and the ring of the native grid: on the same decomposition with the same plans
    if(ck.same_decomposition())
    {
        const vector<double> a = ck.get_rank_blob("sflow_amr_aux",p0->mpirank);
        size_t k = 0;
        auto nx = [&]() { return k<a.size() ? a[k++] : -1.0; };
        bool fits = (int(nx())==(int)gserve.size());
        vector<vector<double>> sv;
        for(size_t l=0; fits && l<gserve.size(); ++l)
        {
            const int n = int(nx());
            fits = (n==(int)gserve[l].size()) && k+n<=a.size();
            if(fits)
            sv.emplace_back(a.begin()+k,a.begin()+k+n);
            if(fits)
            k += n;
        }
        const int nn = fits ? int(nx()) : -1;
        fits = fits && nn==(int)nfill.size() && k+nn<=a.size();
        if(fits)
        {
            for(size_t l=0; l<gserve.size(); ++l)
            for(size_t m=0; m<gserve[l].size(); ++m)
            gserve[l][m].aux = sv[l][m];
            for(int m=0; m<nn; ++m)
            nfill[m].aux = a[k+m];
        }
    }

    return true;
}
