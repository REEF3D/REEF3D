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

#include"fem_solid.h"
#ifdef _OPENMP
#include<omp.h>
#endif
#include<cmath>
#include<cstdint>
#include<algorithm>

// Penalty contact. The normal stiffness is scaled with the nodal mass and the
// element frequency (c_p/h)^2, i.e. it is of the order of the element stiffness
// and stable with the element time step. Friction is regularised Coulomb.

void fem_solid::contact_surface(int i,double pen,const Vec3& n,double mu)
{
    // penalty along the outward normal n of the wall; Coulomb friction with a
    // tangential spring (sticks below mu Fn, slides above)
    // stiffness of the node's own material (rigid bodies: rigid_contact_speed)
    const double c = cnode.empty() ? cp_contact : cnode[i];
    const double w2 = kground*(c/hmin())*(c/hmin());

    const double kn = w2*m[i];
    const double cn = 2.0*contact_zeta*std::sqrt(kn*m[i]);
    const double vn = v[i].dot(n);
    const double Fn = std::max(0.0, kn*pen - cn*vn);

    fcon[i] += Fn*n;

    Vec3 vt = v[i] - vn*n;
    if(plane_strain) vt(1) = 0.0;

    Vec3& u = tspring[i];
    u -= u.dot(n)*n;                    // keep the spring in the tangent plane
    u += dts_cur*vt;
    if(plane_strain) u(1) = 0.0;

    const double kt = kn;
    const double ct = cn;
    Vec3 Ft = -kt*u - ct*vt;
    const double Fmax = mu*Fn;
    const double fn = Ft.norm();
    if(fn>Fmax)
    {
        // sliding: the spring is limited to the Coulomb force
        if(fn>0.0) Ft *= Fmax/fn;
        u = -(Ft + ct*vt)/kt;
        const double un = u.norm(), umax = Fmax/kt;
        if(un>umax && un>0.0) u *= umax/un;
    }
    fcon[i] += Ft;
    touched[i] = 1;
}

void fem_solid::set_bed_sample(int i,double phi,const Vec3& n)
{
    if(bed_ok.size()!=(size_t)nnode())
    clear_bed_samples();
    bed_ok[i] = 1;
    bed_phi[i] = phi;
    bed_n[i] = n;
    bed_x[i] = x[i];
}

void fem_solid::contact_ground()
{
    const Vec3 ez(0.0,0.0,1.0);
    const int nn = nnode();
    if(tspring.size()!=(size_t)nn)
    tspring.assign(nn,Vec3::Zero());
    touched.assign(nn,0);

    // deformable nodes: penalty per node; nodes of rigid bodies: collected for
    // the contact of the whole body (rigid_contact_forces). On threads every thread
    // collects the pairs of its contiguous block of nodes, appended in the order
    // of the blocks: the list is the same as with one thread
    auto& rp = thr_rp;
    rp.resize(nthr);
    for(auto& r : rp) r.v.clear();

    FEM_OMP(omp parallel num_threads(nthr) if(par(nn,PAR_NODES)))
    {
        std::vector<rpair>& rpl = rp[thread_num()].v;
        auto touch = [&](int i,double pen,const Vec3& n,int code)
        {
            const int r = rnode.empty() ? -1 : rnode[i];
            if(r<0)
            contact_surface(i,pen,n,mu_ground);
            else
            {
                rpl.push_back({i,-1,((long long)r<<24) + code,n,pen,v[i].dot(n)});
                touched[i] = 1;
            }
        };

        FEM_OMP(omp for schedule(static))
        for(int i=0; i<nn; ++i)
        {
            if(m[i]>0.0)
            {
                // ground plane: all nodes
                if(ground_on)
                {
                    const double pen = zground - x[i](2);
                    if(pen>0.0)
                    touch(i,pen,ez,0);
                }

                // domain walls and the bed of the fluid grid: free bodies and debris only
                if((!planes.empty() || bed_on) && free_node(i))
                {
                    for(size_t j=0; j<planes.size(); ++j)
                    {
                        const cplane& pl = planes[j];
                        const double pen = pl.d - pl.n.dot(x[i]);
                        if(pen>0.0)
                        touch(i,pen,pl.n,1+(int)std::min<size_t>(j,998));
                    }

                    if(bed_on && !bed_ok.empty() && bed_ok[i])
                    {
                        const double pen = -(bed_phi[i] + bed_n[i].dot(x[i]-bed_x[i]));
                        if(pen>0.0)
                        touch(i,pen,bed_n[i],1001);
                    }
                }
            }

            // nodes out of contact: release the tangential springs
            if(!touched[i])
            tspring[i].setZero();
        }
    }

    for(const auto& r : rp)
    rpairs.insert(rpairs.end(),r.v.begin(),r.v.end());
}

bool fem_solid::share_alive_element(int a,int b) const
{
    for(int q=node_elem_start[a]; q<node_elem_start[a+1]; ++q)
    {
        const element& el = elems[node_elem[q]];
        if(!el.alive)
        continue;
        for(int c=0; c<8; ++c)
        if(el.n[c]==b)
        return true;
    }
    return false;
}

void fem_solid::contact_nodes()
{
    // threads: every thread works on a contiguous block of the candidates and keeps its
    // results in the order of the candidates; they are joined (and the forces applied)
    // in the order of the blocks, i.e. exactly as with one thread
    const int nn = nnode();

    // candidates: nodes on the surface of intact parts (fewer than 8 intact
    // elements) and debris particles
    std::vector<int> cand;
    {
        auto& ct = thr_i;
        ct.resize(nthr);
        for(auto& c : ct) c.v.clear();
        FEM_OMP(omp parallel num_threads(nthr) if(par(nn,PAR_NODES)))
        {
            std::vector<int>& c = ct[thread_num()].v;
            FEM_OMP(omp for schedule(static))
            for(int i=0; i<nn; ++i)
            if(m[i]>0.0 && nalive[i]<8)
            c.push_back(i);
        }
        for(const auto& c : ct)
        cand.insert(cand.end(),c.v.begin(),c.v.end());
    }
    const int nc = (int)cand.size();

    const double d0 = contact_dist*hmin();
    const double cs = d0;           // hash cell size of the pair order

    const double cpmax = cp_contact;
    const double dts = cfl*(dtcrit_el<1.0e30 ? dtcrit_el : dtcrit);

    // neighbour list with a skin: the pairs closer than d0 + skin that can touch
    // (pairs of one intact element or of one rigid body never do); valid while no
    // candidate has moved more than skin/2 and no element has failed since the build
    const double skin = 0.5*hmin();
    const int neroded = n_eroded();
    bool rebuild = (cand!=nl_cand) || nl_x0.size()!=cand.size() || nl_neroded!=neroded || nl_nrigid!=(int)rbs.size();
    if(!rebuild)
    {
        int moved = 0;
        FEM_OMP(omp parallel for schedule(static) num_threads(nthr) if(par(nc,PAR_NODES)) reduction(max:moved))
        for(int q=0; q<nc; ++q)
        if((x[cand[q]]-nl_x0[q]).squaredNorm() > 0.25*skin*skin)
        moved = 1;
        rebuild = moved;
    }

    if(rebuild)
    {
        const double rl = d0 + skin;
        const int64_t B = 1<<20;
        auto key = [&](int ix,int iy,int iz)->int64_t
        {
            return ((int64_t)(ix+B/2)) + B*(((int64_t)(iy+B/2)) + B*((int64_t)(iz+B/2)));
        };

        std::vector<std::pair<int64_t,int>> cells(nc);
        std::vector<int> ci(3*(size_t)nc);
        FEM_OMP(omp parallel for schedule(static) num_threads(nthr) if(par(nc,PAR_NODES)))
        for(int q=0; q<nc; ++q)
        {
            const Vec3& p = x[cand[q]];
            const int ix = (int)std::floor(p(0)/rl), iy = (int)std::floor(p(1)/rl), iz = (int)std::floor(p(2)/rl);
            ci[3*q]=ix; ci[3*q+1]=iy; ci[3*q+2]=iz;
            cells[q] = std::make_pair(key(ix,iy,iz),q);
        }
        std::sort(cells.begin(),cells.end());

        nl_cand = cand;
        nl_x0.resize(nc);
        nl_start.assign(nc+1,0);
        nl_nb.clear();
        auto& nbt = thr_i;
        nbt.resize(nthr);
        for(auto& c : nbt) c.v.clear();
        FEM_OMP(omp parallel num_threads(nthr) if(par(nc,PAR_ELEMS)))
        {
            std::vector<int>& nb = nbt[thread_num()].v;
            FEM_OMP(omp for schedule(static))
            for(int q=0; q<nc; ++q)
            {
                nl_x0[q] = x[cand[q]];
                const int a = cand[q];
                const size_t n0 = nb.size();
                for(int dz=-1; dz<=1; ++dz)
                for(int dy=-1; dy<=1; ++dy)
                for(int dx=-1; dx<=1; ++dx)
                {
                    const int64_t k = key(ci[3*q]+dx,ci[3*q+1]+dy,ci[3*q+2]+dz);
                    auto it = std::lower_bound(cells.begin(),cells.end(),std::make_pair(k,-1));
                    for(; it!=cells.end() && it->first==k; ++it)
                    {
                        const int b = cand[it->second];
                        if(b<=a || (x[a]-x[b]).squaredNorm() >= rl*rl)
                        continue;
                        if(!rnode.empty() && rnode[a]>=0 && rnode[a]==rnode[b])
                        continue;
                        if(share_alive_element(a,b))
                        continue;
                        nb.push_back(b);
                    }
                }
                nl_start[q+1] = (int)(nb.size()-n0);       // count, summed below
            }
        }
        for(int q=0; q<nc; ++q)
        nl_start[q+1] += nl_start[q];
        nl_nb.reserve(nl_start[nc]);
        for(const auto& nb : nbt)
        nl_nb.insert(nl_nb.end(),nb.v.begin(),nb.v.end());
        nl_neroded = neroded;
        nl_nrigid = (int)rbs.size();
        ++nl_builds;
    }

    // the pairs closer than d0, in the order of a hash grid of cell size d0 (the
    // 27 cells around a, then the node number), as without the list
    auto& evt = thr_ev;
    evt.resize(nthr);
    for(auto& e : evt) e.v.clear();
    thr_near.resize(nthr);

    auto contact_pair = [&](std::vector<cevent>& ev,int a,int b,double d)
    {
        const Vec3 r = x[a]-x[b];

        if(share_alive_element(a,b))
        return;

        // nodes of rigid bodies: contact of the whole body
        const int ra = rnode.empty() ? -1 : rnode[a], rb = rnode.empty() ? -1 : rnode[b];
        if(ra>=0 || rb>=0)
        {
            if(ra==rb)
            return;
            // node of the (lower) rigid body first, normal from the partner to it
            int p = a, o = b, rp = ra, ro = rb;
            if(rp<0 || (ro>=0 && ro<rp)) {std::swap(p,o); std::swap(rp,ro);}
            const Vec3 n = (x[p]-x[o])/d;
            const int code = ro>=0 ? 2000+ro : 1002;
            ev.push_back({p,o,Vec3::Zero(),true,{p,o,((long long)rp<<24) + code,n,d0-d,(v[p]-v[o]).dot(n)}});
            return;
        }

        const Vec3 n = r/d;
        const double meff = m[a]*m[b]/(m[a]+m[b]);
        // the softer of the two materials (springs in series)
        const double c = cnode.empty() ? cpmax : std::min(cnode[a],cnode[b]);
        const double kn = kcontact*(c/hmin())*(c/hmin())*meff;
        const double cn = 2.0*contact_zeta*std::sqrt(kn*meff);
        const Vec3 vr = v[a]-v[b];
        const double vn = vr.dot(n);
        const double Fn = std::max(0.0, kn*(d0-d) - cn*vn);

        if(Fn<=0.0)
        return;

        Vec3 F = Fn*n;

        const Vec3 vt = vr - vn*n;
        const double vtn = vt.norm();
        if(vtn>1.0e-12)
        {
            const double ct = 0.5*meff/dts;
            F -= std::min(mu_contact*Fn, ct*vtn)*vt/vtn;
        }

        ev.push_back({a,b,F,false,rpair()});
    };

    FEM_OMP(omp parallel num_threads(nthr) if(par(nc,PAR_ELEMS)))
    {
        std::vector<cevent>& ev = evt[thread_num()].v;
        std::vector<std::pair<int,int>>& near = thr_near[thread_num()].v;   // (cell offset, partner)
        FEM_OMP(omp for schedule(static))
        for(int q=0; q<nc; ++q)
        {
            const int a = cand[q];
            near.clear();
            const Vec3& pa = x[a];
            const int ax = (int)std::floor(pa(0)/cs), ay = (int)std::floor(pa(1)/cs), az = (int)std::floor(pa(2)/cs);
            for(int p=nl_start[q]; p<nl_start[q+1]; ++p)
            {
                const int b = nl_nb[p];
                const Vec3& pb = x[b];
                if((pa-pb).squaredNorm() > 1.000001*d0*d0)     // the exact test follows
                continue;
                const int ox = (int)std::floor(pb(0)/cs)-ax, oy = (int)std::floor(pb(1)/cs)-ay, oz = (int)std::floor(pb(2)/cs)-az;
                if(ox<-1 || ox>1 || oy<-1 || oy>1 || oz<-1 || oz>1)
                continue;
                near.push_back(std::make_pair((oz+1)*9+(oy+1)*3+(ox+1),b));
            }
            std::sort(near.begin(),near.end());
            for(const std::pair<int,int>& c : near)
            {
                const int b = c.second;
                const double d = (x[a]-x[b]).norm();
                if(d>=d0 || d<1.0e-14)
                continue;
                contact_pair(ev,a,b,d);
            }
        }
    }

    // in the order of the candidates
    for(const auto& ev : evt)
    for(const cevent& c : ev.v)
    {
        if(c.rigid)
        rpairs.push_back(c.rp);
        else
        {
            fcon[c.a] += c.F;
            fcon[c.b] -= c.F;
        }
    }
}

void fem_solid::rigid_contact_forces()
{
    // Each rigid body acts as one spring of its effective stiffness k (the
    // debris) against each partner: a wall plane, the ground, the bed, another
    // rigid body (springs in series) or the deformable parts (their own
    // flexibility comes from the FEM). The force of a group is k times the
    // largest penetration, F = k delta_max (+ damping), capped at the crushing
    // force with a permanent set; it is shared by the contact points in
    // proportion to their penetration. A face-on impact at speed u gives the
    // peak u sqrt(k M) after (pi/2) sqrt(M/k).
    if(rpairs.empty())
    {
        rset.clear();
        return;
    }

    std::stable_sort(rpairs.begin(),rpairs.end(),[](const rpair& p,const rpair& q){return p.key<q.key;});
    std::map<std::pair<long long,long long>,double> rnew;
    const long long nn1 = nnode()+1;
    std::vector<double> de;

    size_t q0 = 0;
    while(q0<rpairs.size())
    {
        const long long key = rpairs[q0].key;
        size_t q1 = q0;
        while(q1<rpairs.size() && rpairs[q1].key==key) ++q1;

        const int r = (int)(key>>24), code = (int)(key & 0xFFFFFF);
        const rigid_body& A = rbs[r];
        const bool plane = code<=1001;

        double k = A.k, meff = A.M, cap = A.Fcap, mu = plane ? mu_ground : mu_contact;
        if(code>=2000)
        {
            const rigid_body& B = rbs[code-2000];
            k = A.k*B.k/(A.k+B.k);
            meff = A.M*B.M/(A.M+B.M);
            cap = (A.Fcap>0.0 && B.Fcap>0.0) ? std::min(A.Fcap,B.Fcap) : std::max(A.Fcap,B.Fcap);
        }
        else if(code==1002)
        {
            // the surface of the deformable part acts in series with the debris:
            // penalty of its nodes, of the order of the element stiffness and
            // stable with the element time step (as the contact of deformable parts)
            std::vector<int> part;
            for(size_t q=q0; q<q1; ++q) part.push_back(rpairs[q].b);
            std::sort(part.begin(),part.end());
            part.erase(std::unique(part.begin(),part.end()),part.end());
            double mP = 0.0, KP = 0.0;
            for(int b : part)
            {
                mP += m[b];
                const double c = cnode.empty() ? cp_contact : cnode[b];
                KP += kcontact*(c/hmin())*(c/hmin())*m[b];
            }
            meff = A.M*mP/(A.M+mP);
            if(KP>0.0)
            k = A.k*KP/(A.k+KP);
        }

        // elastic-perfectly plastic spring: crushing beyond the cap leaves a
        // permanent set at the contact points that carry it
        auto pid = [&](const rpair& c){return std::make_pair(key,(long long)c.a*nn1 + (long long)(c.b+1));};
        de.assign(q1-q0,0.0);
        std::vector<double> dp(q1-q0,0.0);
        double S = 0.0, dmax = 0.0;
        for(size_t q=q0; q<q1; ++q)
        {
            auto it = rset.find(pid(rpairs[q]));
            if(it!=rset.end()) dp[q-q0] = it->second;
            de[q-q0] = std::max(0.0, rpairs[q].pen - dp[q-q0]);
            dmax = std::max(dmax,de[q-q0]);
        }
        double Fel = k*dmax;
        if(cap>0.0 && Fel>cap)
        {
            const double dc = cap/k;
            for(size_t q=q0; q<q1; ++q)
            if(de[q-q0]>dc)
            {
                dp[q-q0] = rpairs[q].pen - dc;
                de[q-q0] = dc;
            }
            Fel = cap;
        }
        for(size_t q=q0; q<q1; ++q)
        {
            S += de[q-q0];
            if(dp[q-q0]>0.0) rnew[pid(rpairs[q])] = dp[q-q0];
        }
        if(S<=0.0)
        {
            q0 = q1;
            continue;
        }
        double vG = 0.0;
        for(size_t q=q0; q<q1; ++q)
        vG += de[q-q0]/S*rpairs[q].vn;

        // damping only while the contact unloads: the peak force stays the
        // elastic u sqrt(k M), the rebound loses energy (restitution 0.55 at 50 %)
        const double C = 2.0*debris_zeta*std::sqrt(k*meff);
        double F = std::max(0.0, Fel - C*std::max(vG,0.0));
        if(cap>0.0) F = std::min(F,cap);

        for(size_t q=q0; q<q1; ++q)
        {
            const rpair& c = rpairs[q];
            const double w = de[q-q0]/S;
            const double fn = F*w;
            fcon[c.a] += fn*c.n;
            if(c.b>=0) fcon[c.b] -= fn*c.n;

            if(fn<=0.0)
            continue;

            if(plane)
            {
                // Coulomb friction with a tangential spring (share w of the body stiffness)
                Vec3 vt = v[c.a] - v[c.a].dot(c.n)*c.n;
                if(plane_strain) vt(1) = 0.0;
                Vec3& u = tspring[c.a];
                u -= u.dot(c.n)*c.n;
                u += dts_cur*vt;
                if(plane_strain) u(1) = 0.0;
                const double kt = k*w, ct = C*w;
                Vec3 Ft = -kt*u - ct*vt;
                const double Fmax = mu*fn, ftn = Ft.norm();
                if(ftn>Fmax)
                {
                    if(ftn>0.0) Ft *= Fmax/ftn;
                    u = -(Ft + ct*vt)/kt;
                    const double un = u.norm(), umax = Fmax/kt;
                    if(un>umax && un>0.0) u *= umax/un;
                }
                fcon[c.a] += Ft;
            }
            else
            {
                const Vec3 vr = v[c.a]-v[c.b];
                const Vec3 vt = vr - vr.dot(c.n)*c.n;
                const double vtn = vt.norm();
                if(vtn>1.0e-12)
                {
                    const double mn = m[c.a]*m[c.b]/(m[c.a]+m[c.b]);
                    const Vec3 Ft = -std::min(mu*fn, 0.5*mn/dts_cur*vtn)*vt/vtn;
                    fcon[c.a] += Ft;
                    fcon[c.b] -= Ft;
                }
            }
        }
        q0 = q1;
    }
    rset.swap(rnew);
}
