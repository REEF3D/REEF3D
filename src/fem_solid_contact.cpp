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
    if(tspring.size()!=(size_t)nnode())
    tspring.assign(nnode(),Vec3::Zero());
    touched.assign(nnode(),0);

    // deformable nodes: penalty per node; nodes of rigid bodies: collected for
    // the contact of the whole body (rigid_contact_forces)
    auto touch = [&](int i,double pen,const Vec3& n,int code)
    {
        const int r = rnode.empty() ? -1 : rnode[i];
        if(r<0)
        contact_surface(i,pen,n,mu_ground);
        else
        {
            rpairs.push_back({i,-1,((long long)r<<24) + code,n,pen,v[i].dot(n)});
            touched[i] = 1;
        }
    };

    for(int i=0; i<nnode(); ++i)
    {
        if(m[i]<=0.0)
        continue;

        // ground plane: all nodes
        if(ground_on)
        {
            const double pen = zground - x[i](2);
            if(pen>0.0)
            touch(i,pen,ez,0);
        }

        // domain walls and the bed of the fluid grid: free bodies and debris only
        if(!planes.empty() || bed_on)
        {
            if(!free_node(i))
            continue;

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
    for(int i=0; i<nnode(); ++i)
    if(!touched[i])
    tspring[i].setZero();
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
    // candidates: nodes on the surface of intact parts (fewer than 8 intact
    // elements) and debris particles
    std::vector<int> cand;
    for(int i=0; i<nnode(); ++i)
    if(m[i]>0.0 && nalive[i]<8)
    cand.push_back(i);

    const double d0 = contact_dist*hmin();
    const double cs = d0;           // hash cell size

    const double cpmax = cp_contact;
    const double dts = cfl*(dtcrit_el<1.0e30 ? dtcrit_el : dtcrit);

    const int64_t B = 1<<20;
    auto key = [&](int ix,int iy,int iz)->int64_t
    {
        return ((int64_t)(ix+B/2)) + B*(((int64_t)(iy+B/2)) + B*((int64_t)(iz+B/2)));
    };

    std::vector<std::pair<int64_t,int>> cells(cand.size());
    std::vector<int> ci(3*cand.size());
    for(size_t q=0; q<cand.size(); ++q)
    {
        const Vec3& p = x[cand[q]];
        const int ix = (int)std::floor(p(0)/cs), iy = (int)std::floor(p(1)/cs), iz = (int)std::floor(p(2)/cs);
        ci[3*q]=ix; ci[3*q+1]=iy; ci[3*q+2]=iz;
        cells[q] = std::make_pair(key(ix,iy,iz),cand[q]);
    }
    std::sort(cells.begin(),cells.end());

    for(size_t q=0; q<cand.size(); ++q)
    {
        const int a = cand[q];

        for(int dz=-1; dz<=1; ++dz)
        for(int dy=-1; dy<=1; ++dy)
        for(int dx=-1; dx<=1; ++dx)
        {
            const int64_t k = key(ci[3*q]+dx,ci[3*q+1]+dy,ci[3*q+2]+dz);
            auto it = std::lower_bound(cells.begin(),cells.end(),std::make_pair(k,-1));

            for(; it!=cells.end() && it->first==k; ++it)
            {
                const int b = it->second;
                if(b<=a)
                continue;

                Vec3 r = x[a]-x[b];
                const double d = r.norm();
                if(d>=d0 || d<1.0e-14)
                continue;

                if(share_alive_element(a,b))
                continue;

                // nodes of rigid bodies: contact of the whole body
                const int ra = rnode.empty() ? -1 : rnode[a], rb = rnode.empty() ? -1 : rnode[b];
                if(ra>=0 || rb>=0)
                {
                    if(ra==rb)
                    continue;
                    // node of the (lower) rigid body first, normal from the partner to it
                    int p = a, o = b, rp = ra, ro = rb;
                    if(rp<0 || (ro>=0 && ro<rp)) {std::swap(p,o); std::swap(rp,ro);}
                    const Vec3 n = (x[p]-x[o])/d;
                    const int code = ro>=0 ? 2000+ro : 1002;
                    rpairs.push_back({p,o,((long long)rp<<24) + code,n,d0-d,(v[p]-v[o]).dot(n)});
                    continue;
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
                continue;

                Vec3 F = Fn*n;

                const Vec3 vt = vr - vn*n;
                const double vtn = vt.norm();
                if(vtn>1.0e-12)
                {
                    const double ct = 0.5*meff/dts;
                    F -= std::min(mu_contact*Fn, ct*vtn)*vt/vtn;
                }

                fcon[a] += F;
                fcon[b] -= F;
            }
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
