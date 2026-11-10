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

#include"cfd_pocket.h"
#include"cfd_singlephase.h"
#include"cfd_amr.h"
#include"momentum_rk.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"increment.h"
#include"reef3d_mpmd.h"
#include<mpi.h>
#include<map>
#include<set>
#include<utility>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<cstdlib>

namespace
{
// a pocket of the last labelling: mass (its volume at p_ref), volume, pressure, volume of its level set
struct pocket
{
    double m, Vk, p, V;
};

// the pockets of the last labelling (the same on all ranks), the step that followed it
std::map<double,pocket> last;
double stamp = 0.0;
double dt_last = 0.0;
std::ofstream plog;

// the labels of the cells: >= 0 a cell in the air (its component), WATER a leaf cell in the water,
// -1 none (not a leaf cell of the grid: covered, solid, outside)
const double WATER = -3.0;

inline int fid(lexer *p, int i, int j, int k)
{
    return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
}

struct dsu
{
    std::vector<int> u;
    explicit dsu(int n) : u(n) { for(int a=0; a<n; ++a) u[a] = a; }
    int find(int a) { while(u[a]!=a) { u[a] = u[u[a]]; a = u[a]; } return a; }
    void join(int a, int b) { a = find(a); b = find(b); if(a!=b) u[std::max(a,b)] = std::min(a,b); }
};

// the widths of a cell
inline void cell_w(lexer *p, int i, int j, int k, double *w)
{
    const int m = increment::marge;
    w[0] = p->DXN[i+m];
    w[1] = p->DYN[j+m];
    w[2] = p->DZN[k+m];
}

// the stats of a component: volume, that of its cells in the air, their number, the flux of the
// water out of it (dV/dt) and the area of its faces to the water, the overlaps of its cells in the
// air with the pockets of the last labelling (-1: the atmosphere)
struct comp
{
    double V = 0.0, Va = 0.0, n = 0.0, Q = 0.0, A = 0.0;
    std::map<double,double> ov;
};
}

void cfd_pocket_update(lexer *p, fdm *a, ghostcell *pgc, momentum *pmom, cfd_amr *pamr)
{
    if(p->F41!=2)
    return;

    if(pamr!=nullptr)
    {
        pamr->pocket_update(p,pgc);
        return;
    }

    momentum_rk *m = dynamic_cast<momentum_rk*>(pmom);
    if(m==nullptr || m->single()==nullptr)
    return;
    cfd_singlephase *s = m->single();

    std::vector<cfd_pocket_grid> G(1);
    G[0].p = p;
    G[0].a = a;
    G[0].s = s;
    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    if(p->flag4[fid(p,i,j,k)]>0)
    G[0].cells.push_back({i,j,k});

    const std::vector<cfd_pocket_edge> E;
    const std::function<void(int)> ex = [&](int kind)
    {
        if(kind==0)
        pgc->gcparax(p,s->lab,4);
        else
        {
            pgc->gcparax(p,s->pid,4);
            pgc->gcparax(p,s->pa,4);
        }
    };

    cfd_pocket_core(p,G,E,ex);
}

void cfd_pocket_core(lexer *p0, std::vector<cfd_pocket_grid> &G, const std::vector<cfd_pocket_edge> &E, const std::function<void(int)> &exchange)
{
    MPI_Comm comm = reef3d_comm();
    int rank, size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
    const int ng = (int)G.size();
    const int jd = p0->j_dir==1 ? 1 : 0;
    const int nb[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

    // not the cells inside the solids and the topography of the direct forcing (G 505 0): they are
    // cells of the grid, their level set extended, but no air of a pocket (the deck above the bays
    // would join them to the atmosphere)
    for(int g=0; g<ng; ++g)
    {
        lexer *p = G[g].p;
        fdm *a = G[g].a;
        if(p->G505!=0 || (p->solidread<=0 && p->topoforcing<=0))
        continue;
        auto &cl = G[g].cells;
        cl.erase(std::remove_if(cl.begin(),cl.end(),[&](const std::array<int,3> &c)
                 {
                     return (p->solidread>0 && a->solid(c[0],c[1],c[2])<0.0) || (p->topoforcing>0 && a->topo(c[0],c[1],c[2])<0.0);
                 }),cl.end());
    }

    // the cells in the air of this rank, joined within the grids and across the patch edges
    std::vector<std::vector<int>> node(ng);
    int nn = 0;
    for(int g=0; g<ng; ++g)
    {
        lexer *p = G[g].p;
        node[g].assign((size_t)p->imax*p->jmax*p->kmax,-1);
        for(const auto &c : G[g].cells)
        if(G[g].a->phi(c[0],c[1],c[2])<0.0)
        node[g][fid(p,c[0],c[1],c[2])] = nn++;
    }

    dsu D(nn);
    for(int g=0; g<ng; ++g)
    {
        lexer *p = G[g].p;
        for(const auto &c : G[g].cells)
        {
            const int n = node[g][fid(p,c[0],c[1],c[2])];
            if(n<0)
            continue;
            for(int s=0; s<6; s+=2)
            {
                if(s==2 && jd==0)
                continue;
                const int m = node[g][fid(p,c[0]+nb[s][0],c[1]+nb[s][1],c[2]+nb[s][2])];
                if(m>=0)
                D.join(n,m);
            }
        }
    }
    for(const cfd_pocket_edge &e : E)
    {
        const int n1 = node[e.g1][fid(G[e.g1].p,e.c1[0],e.c1[1],e.c1[2])];
        const int n2 = node[e.g2][fid(G[e.g2].p,e.c2[0],e.c2[1],e.c2[2])];
        if(n1>=0 && n2>=0)
        D.join(n1,n2);
    }

    // labels unique over the ranks
    std::vector<int> cidx(nn,-1);
    long long nc = 0;
    for(int n=0; n<nn; ++n)
    if(D.find(n)==n)
    cidx[n] = (int)nc++;
    long long off = 0;
    MPI_Exscan(&nc,&off,1,MPI_LONG_LONG,MPI_SUM,comm);
    if(rank==0)
    off = 0;

    for(int g=0; g<ng; ++g)
    {
        lexer *p = G[g].p;
        field &L = G[g].s->lab;
        const int sz = p->imax*p->jmax*p->kmax;
        for(int q=0; q<sz; ++q)
        L.V[q] = -1.0;
        for(const auto &c : G[g].cells)
        {
            const int n = node[g][fid(p,c[0],c[1],c[2])];
            L(c[0],c[1],c[2]) = (n>=0) ? double(off + cidx[D.find(n)]) : WATER;
        }
    }

    exchange(0);

    // the labels of the neighbours on the other ranks and the sibling patches: pairs of equal ones
    std::set<std::pair<long long,long long>> pairs;
    for(int g=0; g<ng; ++g)
    {
        field &L = G[g].s->lab;
        for(const auto &c : G[g].cells)
        {
            const double l0 = L(c[0],c[1],c[2]);
            if(l0<0.0)
            continue;
            for(int s=0; s<6; ++s)
            {
                if((s==2 || s==3) && jd==0)
                continue;
                const double l1 = L(c[0]+nb[s][0],c[1]+nb[s][1],c[2]+nb[s][2]);
                if(l1>=0.0 && l1!=l0)
                pairs.insert(std::make_pair((long long)std::min(l0,l1),(long long)std::max(l0,l1)));
            }
        }
    }

    std::vector<long long> sp;
    for(const auto &q : pairs)
    {
        sp.push_back(q.first);
        sp.push_back(q.second);
    }
    int ns = (int)sp.size();
    std::vector<int> cnt(size), dsp(size,0);
    MPI_Allgather(&ns,1,MPI_INT,cnt.data(),1,MPI_INT,comm);
    for(int r=1; r<size; ++r)
    dsp[r] = dsp[r-1]+cnt[r-1];
    std::vector<long long> ap(dsp[size-1]+cnt[size-1]);
    MPI_Allgatherv(sp.data(),ns,MPI_LONG_LONG,ap.data(),cnt.data(),dsp.data(),MPI_LONG_LONG,comm);

    // the components over all ranks: the smallest of the connected labels
    std::map<long long,long long> up;
    auto rep = [&](long long l) -> long long
    {
        long long r = l;
        for(auto it=up.find(r); it!=up.end() && it->second!=r; it=up.find(r))
        r = it->second;
        for(auto it=up.find(l); it!=up.end() && it->second!=r; it=up.find(l))
        {
            l = it->second;
            it->second = r;
        }
        return r;
    };
    for(size_t q=0; q+1<ap.size(); q+=2)
    {
        const long long a = rep(ap[q]), b = rep(ap[q+1]);
        if(a!=b)
        up[std::max(a,b)] = std::min(a,b);
    }

    // per component: the volumes, the overlaps with the pockets of the last labelling, the flux of
    // the water across its faces to the cells in the water
    std::map<long long,comp> S;
    for(int g=0; g<ng; ++g)
    {
        lexer *p = G[g].p;
        fdm *a = G[g].a;
        field &L = G[g].s->lab;
        field &pid = G[g].s->pid;
        field *vel[3] = {&a->u,&a->v,&a->w};
        for(const auto &c : G[g].cells)
        {
            const int i = c[0], j = c[1], k = c[2];
            const double ph = a->phi(i,j,k);
            double w[3];
            cell_w(p,i,j,k,w);
            const double h = (jd==1) ? (w[0]+w[1]+w[2])/3.0 : 0.5*(w[0]+w[2]);
            const double dv = w[0]*w[1]*w[2];

            if(ph<0.0)
            {
                comp &C = S[rep((long long)L(i,j,k))];
                const double fr = std::min(1.0,0.5-ph/h);
                C.V += fr*dv;
                C.Va += fr*dv;
                C.n += 1.0;
                const double o = pid(i,j,k);
                if(o>=-1.0)
                C.ov[o] += fr*dv;

                // the faces to the cells in the water: u.n out of the pocket
                for(int s=0; s<6; ++s)
                {
                    const int d = s/2;
                    if(d==1 && jd==0)
                    continue;
                    if(L(i+nb[s][0],j+nb[s][1],k+nb[s][2])!=WATER)
                    continue;
                    const double af = dv/w[d];
                    const double un = (s%2==0) ? (*vel[d])(i,j,k) : -(*vel[d])(i+nb[s][0],j+nb[s][1],k+nb[s][2]);
                    C.Q += un*af;
                    C.A += af;
                }
                continue;
            }

            // a cell in the water at the surface: its fraction to the component of a neighbour
            if(ph>=0.5*h)
            continue;
            for(int s=0; s<6; ++s)
            {
                if((s==2 || s==3) && jd==0)
                continue;
                const double l1 = L(i+nb[s][0],j+nb[s][1],k+nb[s][2]);
                if(l1>=0.0)
                {
                    S[rep((long long)l1)].V += (0.5-ph/h)*dv;
                    break;
                }
            }
        }
    }

    // the faces at the patch edges between a cell in the air and one in the water
    for(const cfd_pocket_edge &e : E)
    {
        const double l1 = G[e.g1].s->lab(e.c1[0],e.c1[1],e.c1[2]);
        const double l2 = G[e.g2].s->lab(e.c2[0],e.c2[1],e.c2[2]);
        if(l1>=0.0 && l2==WATER)
        {
            comp &C = S[rep((long long)l1)];
            C.Q += double(e.dir)*e.uf*e.af;
            C.A += e.af;
        }
        if(l2>=0.0 && l1==WATER)
        {
            comp &C = S[rep((long long)l2)];
            C.Q -= double(e.dir)*e.uf*e.af;
            C.A += e.af;
        }
    }

    // all components to all ranks
    std::vector<double> sd;
    for(const auto &kv : S)
    {
        sd.push_back(double(kv.first));
        sd.push_back(kv.second.V);
        sd.push_back(kv.second.Va);
        sd.push_back(kv.second.n);
        sd.push_back(kv.second.Q);
        sd.push_back(kv.second.A);
        sd.push_back(double(kv.second.ov.size()));
        for(const auto &o : kv.second.ov)
        {
            sd.push_back(o.first);
            sd.push_back(o.second);
        }
    }
    int nd = (int)sd.size();
    MPI_Allgather(&nd,1,MPI_INT,cnt.data(),1,MPI_INT,comm);
    dsp[0] = 0;
    for(int r=1; r<size; ++r)
    dsp[r] = dsp[r-1]+cnt[r-1];
    std::vector<double> ad(dsp[size-1]+cnt[size-1]);
    MPI_Allgatherv(sd.data(),nd,MPI_DOUBLE,ad.data(),cnt.data(),dsp.data(),MPI_DOUBLE,comm);

    std::map<long long,comp> T;
    for(size_t q=0; q<ad.size();)
    {
        comp &C = T[(long long)ad[q]];
        C.V += ad[q+1];
        C.Va += ad[q+2];
        C.n += ad[q+3];
        C.Q += ad[q+4];
        C.A += ad[q+5];
        const int no = (int)ad[q+6];
        q += 7;
        for(int o=0; o<no; ++o)
        {
            C.ov[ad[q]] += ad[q+1];
            q += 2;
        }
    }

    // the atmosphere: the largest component and the small ones
    long long atm = -1;
    double vmax = -1.0;
    for(const auto &kv : T)
    if(kv.second.V>vmax)
    {
        vmax = kv.second.V;
        atm = kv.first;
    }
    const double nmin = (jd==1) ? 125.0 : 25.0;

    // the overlaps of each pocket of the last labelling over all components
    std::map<double,double> tot;
    for(const auto &kv : T)
    for(const auto &o : kv.second.ov)
    if(o.first>=0.0)
    tot[o.first] += o.second;

    // the pockets: their mass, volume and pressure
    const double pref = p0->W36, gam = p0->W33, row = p0->W1, dt = p0->dt;
    stamp += 1.0;
    std::map<long long,std::pair<double,double>> pk;     // component -> id, pressure
    std::map<double,pocket> now;
    double idx = 0.0;
    for(const auto &kv : T)
    {
        const comp &C = kv.second;
        if(kv.first==atm || C.n<nmin || C.V<=0.0)
        continue;

        // the shares of the pockets of the last labelling; the air of the atmosphere (or of a pocket
        // no longer known) at p_ref, its volume scaled to the whole volume (with the cells in the
        // water at the surface). The cells that were in the water bring no air
        double m = 0.0, Vk = 0.0, mp = 0.0;
        bool old = false;
        for(const auto &o : C.ov)
        {
            const auto it = (o.first>=0.0) ? last.find(o.first) : last.end();
            if(it!=last.end() && tot[o.first]>0.0)
            {
                const double s = o.second/tot[o.first];
                m += s*it->second.m;
                Vk += s*it->second.Vk;
                mp += s*it->second.m*it->second.p;
                old = true;
            }
            else
            {
                const double v = o.second*C.V/std::max(C.Va,1.0e-30);
                m += v;
                Vk += v;
            }
        }

        double pr = 0.0;
        if(!old || m<=0.0)
        {
            // a new pocket: the volume of its level set at p_ref
            m = Vk = C.V;
        }
        else
        {
            // the volume over the last step from the flux of the water, drawn towards the volume of the
            // level set by 2 % per step (the flux and the level set move the surface a little
            // differently, a drift of some % within a second on 8 cells across); bounded against the
            // level set (a factor 4: far beyond the compression of an air pocket in water waves)
            Vk += dt_last*C.Q;
            Vk += 0.02*(C.V-Vk);
            Vk = std::max(0.25*C.V,std::min(4.0*C.V,Vk));

            // the gain of the loop pressure -> flux -> volume -> pressure over a step, estimated with
            // the water around the pocket moved over the length 2 V / A; the pressure follows its
            // target with 1/(1+G0): implicit for a small pocket at a large step (G0 > 1), explicit
            // for a large one (G0 << 1). (A larger factor, 4 G0 or adaptive, lags the pressure behind
            // the volume and drove the oscillations of the pockets under the deck of validation 89)
            const double p1 = mp/m;
            const double pt = pref*(std::pow(m/Vk,gam) - 1.0);
            const double Gl = gam*(p1+pref)*dt*dt*C.A*C.A/(2.0*row*Vk*Vk);
            pr = p1 + (pt-p1)/(1.0+Gl);

            if(rank==0 && std::getenv("REEF3D_POCKET_DIAG")!=nullptr)
            {
                std::cout<<"POCKET t "<<p0->simtime<<" n "<<C.n<<" V "<<C.V<<" Va "<<C.Va<<" Vk "<<Vk<<" dtQ "<<dt_last*C.Q<<" A "<<C.A<<" m "<<m<<" p1 "<<p1<<" pt "<<pt<<" G "<<Gl<<" p "<<pr<<" ov";
                for(const auto &o : C.ov)
                std::cout<<" "<<o.first<<":"<<o.second;
                std::cout<<std::endl;
            }
        }

        const double id = stamp*1.0e6 + idx;
        idx += 1.0;
        pk[kv.first] = std::make_pair(id,pr);
        now[id] = pocket{m,Vk,pr,C.V};
    }
    last.swap(now);
    dt_last = dt;

    // the pocket and its pressure in the cells
    for(int g=0; g<ng; ++g)
    {
        fdm *a = G[g].a;
        field &L = G[g].s->lab;
        field &pid = G[g].s->pid;
        field &pa = G[g].s->pa;
        for(const auto &c : G[g].cells)
        {
            const int i = c[0], j = c[1], k = c[2];
            if(a->phi(i,j,k)>=0.0)
            {
                pid(i,j,k) = -2.0;
                pa(i,j,k) = 0.0;
                continue;
            }
            const auto it = pk.find(rep((long long)L(i,j,k)));
            if(it==pk.end())
            {
                pid(i,j,k) = -1.0;
                pa(i,j,k) = 0.0;
            }
            else
            {
                pid(i,j,k) = it->second.first;
                pa(i,j,k) = it->second.second;
            }
        }
    }

    exchange(1);

    // the log: the number of pockets, the five largest
    if(rank==0)
    {
        if(!plog.is_open())
        {
            plog.open("./REEF3D_Log/REEF3D_CFD_Pockets.dat");
            plog<<"# F 41 2: the closed air pockets. time, number, then for the five largest: volume of the level set [m3], volume [m3], pressure [Pa]"<<std::endl;
        }
        std::vector<pocket> v;
        for(const auto &kv : last)
        v.push_back(kv.second);
        std::sort(v.begin(),v.end(),[](const pocket &x, const pocket &y) { return x.Vk>y.Vk; });
        plog<<std::setprecision(8)<<p0->simtime<<" "<<v.size();
        for(size_t q=0; q<std::min<size_t>(5,v.size()); ++q)
        plog<<" "<<v[q].V<<" "<<v[q].Vk<<" "<<v[q].p;
        plog<<std::endl;

        if(p0->count%p0->P12==0 && !v.empty())
        std::cout<<"F 41 2: "<<v.size()<<" air pockets, the largest "<<std::setprecision(4)<<v[0].Vk<<" m3 at "<<v[0].p<<" Pa"<<std::endl;
    }
}
