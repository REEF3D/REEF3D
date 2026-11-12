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
#include"cfd_pocket.h"
#include"cfd_singlephase.h"
#include"momentum_rk.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// F 41 2: the grids of the labelling of the air pockets: the leaf cells of level 0 and of the patches
// (lq), the faces at the patch edges (cf: the fine cell and the coarse leaf cell), the ghost cells
// of level 0 from the other ranks (gcparax) and those of the patches from their siblings (fill
// with siblings_only: not across the patch edges)
void cfd_amr::pocket_update(lexer *p, ghostcell *pgc)
{
    std::vector<cfd_pocket_grid> G;
    for(int g=-1; g<npat; ++g)
    {
        cfd_pocket_grid x;
        x.p = glex(g);
        x.a = gfd(g);
        x.s = gmom(g)->single();
        if(x.s==nullptr)
        return;
        const int jk = x.p->jmax*x.p->kmax;
        x.cells.reserve(lq[g+1].size());
        for(int qq : lq[g+1])
        {
            const int i = qq/jk + x.p->imin;
            const int r = qq%jk;
            x.cells.push_back({i, r/x.p->kmax + x.p->jmin, r%x.p->kmax + x.p->kmin});
        }
        G.push_back(std::move(x));
    }

    std::vector<cfd_pocket_edge> E;
    E.reserve(cf.size());
    for(const cfd_amr_cf &F : cf)
    {
        cfd_pocket_edge e;
        e.g1 = F.pid+1;
        e.g2 = F.gc+1;
        for(int d=0; d<3; ++d)
        {
            e.c1[d] = F.fi[d];
            e.c2[d] = F.ci[d];
        }
        e.d = F.d;
        e.dir = F.dir;
        cfd_amr_patch *c = CP(F.pid);
        fdm *a = c->a;
        lexer *pp = c->pp;
        field &uf = (F.d==0) ? (field&)a->u : (F.d==1 ? (field&)a->v : (field&)a->w);
        e.uf = uf(F.ff[0],F.ff[1],F.ff[2]);
        const int m = marge;
        const double w[3] = {pp->DXN[F.fi[0]+m],pp->DYN[F.fi[1]+m],pp->DZN[F.fi[2]+m]};
        e.af = 1.0;
        for(int t=0; t<3; ++t)
        if(t!=F.d)
        e.af *= w[t];
        E.push_back(e);
    }

    const std::function<void(int)> ex = [&](int kind)
    {
        std::vector<std::function<field&(int)>> fl;
        if(kind==0)
        fl.push_back([this](int g) -> field& { return gmom(g)->single()->lab; });
        else
        {
            fl.push_back([this](int g) -> field& { return gmom(g)->single()->pid; });
            fl.push_back([this](int g) -> field& { return gmom(g)->single()->pa; });
        }

        for(auto &f : fl)
        pgc->gcparax(p0,f(-1),4);

        vector<fspec> fs;
        for(auto &f : fl)
        {
            fspec s;
            s.type = CELL;
            s.mode = 2;
            s.f = f;
            fs.push_back(s);
        }
        for(int l=1; l<=maxlev; ++l)
        fill(l,fs,true);
    };

    cfd_pocket_core(p0,G,E,ex);
}
