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

void fem_solid::contact_ground()
{
    double cpmax = 0.0;
    for(const material& mt : mats)
    cpmax = std::max(cpmax,mt.cp);
    const double w2 = kground*(cpmax/hmin())*(cpmax/hmin());
    const double dts = cfl*dtcrit;

    for(int i=0; i<nnode(); ++i)
    {
        if(m[i]<=0.0)
        continue;

        const double pen = zground - x[i](2);
        if(pen<=0.0)
        continue;

        const double kn = w2*m[i];
        const double cn = 2.0*contact_zeta*std::sqrt(kn*m[i]);
        const double vn = v[i](2);
        const double Fn = std::max(0.0, kn*pen - cn*vn);

        Vec3 vt = v[i];
        vt(2) = 0.0;
        const double vtn = vt.norm();

        fcon[i](2) += Fn;

        if(vtn>1.0e-12)
        {
            const double ct = 0.5*m[i]/dts;
            const double Ft = std::min(mu_ground*Fn, ct*vtn);
            fcon[i] -= Ft*vt/vtn;
        }
    }
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

    double cpmax = 0.0;
    for(const material& mt : mats)
    cpmax = std::max(cpmax,mt.cp);
    const double w2 = kcontact*(cpmax/hmin())*(cpmax/hmin());
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
                const double kn = w2*meff;
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
