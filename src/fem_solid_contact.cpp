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

    for(int i=0; i<nnode(); ++i)
    {
        if(m[i]<=0.0)
        continue;

        // ground plane: all nodes
        if(ground_on)
        {
            const double pen = zground - x[i](2);
            if(pen>0.0)
            contact_surface(i,pen,ez,mu_ground);
        }

        // domain walls and the bed of the fluid grid: free bodies and debris only
        if(!planes.empty() || bed_on)
        {
            if(!free_node(i))
            continue;

            for(const cplane& pl : planes)
            {
                const double pen = pl.d - pl.n.dot(x[i]);
                if(pen>0.0)
                contact_surface(i,pen,pl.n,mu_ground);
            }

            if(bed_on && !bed_ok.empty() && bed_ok[i])
            {
                const double pen = -(bed_phi[i] + bed_n[i].dot(x[i]-bed_x[i]));
                if(pen>0.0)
                contact_surface(i,pen,bed_n[i],mu_ground);
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
    const double dts = cfl*dtcrit;

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
