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

#include"fem_solid.h"
#include<cmath>
#include<stdexcept>
#include<algorithm>
#include<sstream>

fem_solid::fem_solid()
{
}

fem_solid::~fem_solid()
{
}

void fem_solid::set_lattice(double ox_,double oy_,double oz_,double hx_,double hy_,double hz_)
{
    ox=ox_; oy=oy_; oz=oz_;
    hx=hx_; hy=hy_; hz=hz_;
}

int fem_solid::add_material(const material& mt)
{
    material q = mt;
    q.lambda = q.E*q.nu/((1.0+q.nu)*(1.0-2.0*q.nu));
    q.mu = q.E/(2.0*(1.0+q.nu));
    q.cp = std::sqrt((q.lambda+2.0*q.mu)/q.rho);

    for(size_t n=0; n<mats.size(); ++n)
    if(mats[n].id==q.id)
    {
        mats[n] = q;
        return (int)n;
    }

    mats.push_back(q);
    return (int)mats.size()-1;
}

void fem_solid::add_box(double x0,double x1,double y0,double y1,double z0,double z1,int matid)
{
    shape_box b = {x0,x1,y0,y1,z0,z1,matid,false};
    boxes.push_back(b);
    shapes.push_back({0,(int)boxes.size()-1});
}

void fem_solid::remove_box(double x0,double x1,double y0,double y1,double z0,double z1)
{
    shape_box b = {x0,x1,y0,y1,z0,z1,-1,true};
    boxes.push_back(b);
    shapes.push_back({0,(int)boxes.size()-1});
}

void fem_solid::add_fix(double x0,double x1,double y0,double y1,double z0,double z1,bool fx,bool fy,bool fz)
{
    fix_box f = {x0,x1,y0,y1,z0,z1,{fx,fy,fz}};
    fixes.push_back(f);
}

// ----------------------------------------------------------------------
// setup
// ----------------------------------------------------------------------

void fem_solid::build()
{
    if(mats.empty())
    throw std::runtime_error("FEM: no material defined");

    if(hx<=0.0 || hy<=0.0 || hz<=0.0)
    throw std::runtime_error("FEM: lattice spacing missing or not positive (keyword 'lattice')");

    voxelise();
    make_nodes_elements();

    if(elems.empty())
    throw std::runtime_error("FEM: the geometry does not contain any voxel");

    shape_derivatives();
    compute_mass();

    // fixed dofs
    fixed.assign(nnode(),0);
    for(const fix_box& f : fixes)
    for(int i=0; i<nnode(); ++i)
    {
        const Vec3& p = X[i];
        const double eps = 1.0e-9*hmin();
        if(p(0)>=f.x0-eps && p(0)<=f.x1+eps && p(1)>=f.y0-eps && p(1)<=f.y1+eps && p(2)>=f.z0-eps && p(2)<=f.z1+eps)
        for(int d=0; d<3; ++d)
        if(f.f[d]) fixed[i] |= (unsigned char)(1<<d);
    }

    // monitors: nearest node
    mons.clear();
    for(auto& mp : mon_pts)
    {
        int best = 0;
        double dbest = 1.0e300;
        for(int i=0; i<nnode(); ++i)
        {
            double d = (X[i]-mp.second).squaredNorm();
            if(d<dbest) {dbest = d; best = i;}
        }
        mons.push_back({mp.first,best});
    }

    // critical time step of the explicit scheme (hex8, lumped mass)
    dtcrit = 1.0e30;
    for(const element& e : elems)
    dtcrit = std::min(dtcrit, hmin()/mats[e.mat].cp);

    // crack band check: the softening branch must not snap back
    for(const material& mt : mats)
    if(mt.type==MAT_CONCRETE)
    {
        double e0 = mt.ft/mt.E;
        if(mt.ft>0.0 && mt.Gf/(mt.ft*helem) <= 0.5*e0)
        throw std::runtime_error("FEM: concrete material "+std::to_string(mt.id)+": elements too large for the crack band (h > 2 E Gf / ft^2), refine the lattice or increase Gf");
        double e0c = mt.fc/mt.E;
        if(mt.fc>0.0 && mt.Gc/(mt.fc*helem) <= 0.5*e0c)
        throw std::runtime_error("FEM: concrete material "+std::to_string(mt.id)+": elements too large for the compressive crack band (h > 2 E Gc / fc^2)");
    }

    x = X;
    v.assign(nnode(),Vec3::Zero());
    fint.assign(nnode(),Vec3::Zero());
    fext.assign(nnode(),Vec3::Zero());
    fcon.assign(nnode(),Vec3::Zero());
    Mcpl.assign(nnode(),0.0);
    mfl.assign(nnode(),0.0);
    Mu_cpl.assign(nnode(),Vec3::Zero());
    fcpl.assign(nnode(),Vec3::Zero());

    build_surface();
    count_bodies();

    t = 0.0;
    built = true;
}

void fem_solid::compute_mass()
{
    m.assign(nnode(),0.0);
    vnode.assign(nnode(),0.0);
    Vel = hx*hy*hz;
    helem = std::cbrt(Vel);

    for(const element& e : elems)
    {
        const double me = mats[e.mat].rho*Vel/8.0;
        for(int a=0; a<8; ++a)
        {
            m[e.n[a]] += me;
            vnode[e.n[a]] += Vel/8.0;
        }
    }
}

// ----------------------------------------------------------------------
// time stepping
// ----------------------------------------------------------------------

void fem_solid::clear_loads()
{
    std::fill(fext.begin(),fext.end(),Vec3::Zero());
}

void fem_solid::clear_coupling()
{
    std::fill(Mcpl.begin(),Mcpl.end(),0.0);
    std::fill(mfl.begin(),mfl.end(),0.0);
    std::fill(Mu_cpl.begin(),Mu_cpl.end(),Vec3::Zero());
}

double fem_solid::extent() const
{
    Vec3 lo = Vec3::Constant(1.0e300), hi = Vec3::Constant(-1.0e300);
    for(const Vec3& p : X)
    {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    return X.empty() ? 0.0 : (hi-lo).maxCoeff();
}

void fem_solid::advance(double dt)
{
    if(!built)
    throw std::runtime_error("FEM: advance() before build()");

    if(dt<=0.0)
    return;

    // effective masses: the enclosed fluid reduces the mass that is accelerated
    // by the structural forces (at most to 10 %), the time step follows
    std::vector<double> mr(nnode());
    double rmin = 1.0;
    for(int i=0; i<nnode(); ++i)
    {
        mr[i] = std::max(m[i]-mfl[i], 0.1*m[i]);
        if(m[i]>0.0)
        rmin = std::min(rmin, mr[i]/m[i]);
    }

    const int nsub = std::max(1,(int)std::ceil(dt/(cfl*dtcrit*std::sqrt(rmin)) - 1.0e-9));
    const double dts = dt/double(nsub);
    nsub_last = nsub;

    // merge the attached fluid parcels with the nodes (inelastic), the
    // parcels move with the nodes during the step
    std::vector<double> mt(nnode());
    for(int i=0; i<nnode(); ++i)
    {
        mt[i] = mr[i] + Mcpl[i];
        if(Mcpl[i]>0.0 && mt[i]>0.0)
        v[i] = (mr[i]*v[i] + Mu_cpl[i])/mt[i];

        if(fixed[i])
        for(int d=0; d<3; ++d)
        if(fixed[i] & (1<<d)) v[i](d) = 0.0;
        if(plane_strain)
        v[i](1) = 0.0;

    }

    for(int s=0; s<nsub; ++s)
    {
        std::fill(fint.begin(),fint.end(),Vec3::Zero());
        std::fill(fcon.begin(),fcon.end(),Vec3::Zero());

        internal_forces(dts);

        if(ground_on)
        contact_ground();

        if(contact_on && (nbodies>1 || surf_dirty || n_eroded()>0))
        contact_nodes();

        const double alpha = alpha_damp + (t<relax_time ? relax_alpha : 0.0);

        for(int i=0; i<nnode(); ++i)
        {
            if(m[i]<=0.0)
            continue;

            const Vec3 F = fext[i] + fcon[i] - fint[i] + mr[i]*grav;

            v[i] += dts*(F/mt[i] - alpha*v[i]);

            if(fixed[i])
            for(int d=0; d<3; ++d)
            if(fixed[i] & (1<<d)) v[i](d) = 0.0;

            if(plane_strain)
            v[i](1) = 0.0;

            if(!v[i].allFinite())
            {
                std::ostringstream os;
                os<<"FEM: non-finite velocity at node "<<i<<" x "<<x[i].transpose()<<" m "<<m[i]<<" m_f "<<mfl[i]
                  <<" M_cpl "<<Mcpl[i]<<" M u_f "<<Mu_cpl[i].transpose()<<" f_ext "<<fext[i].transpose()<<" f_int "<<fint[i].transpose()
                  <<" f_con "<<fcon[i].transpose()<<" intact elements "<<nalive[i];
                throw std::runtime_error(os.str());
            }

            x[i] += dts*v[i];
        }

        t += dts;
    }

    // force of the attached fluid on the nodes over the step: the parcels
    // start with u_f and end with the node velocity
    for(int i=0; i<nnode(); ++i)
    fcpl[i] = (Mcpl[i]>0.0) ? Vec3((Mu_cpl[i] - Mcpl[i]*v[i])/dt) : Vec3::Zero();

    if(surf_dirty)
    {
        build_surface();
        count_bodies();
    }
}

// ----------------------------------------------------------------------
// diagnostics
// ----------------------------------------------------------------------

int fem_solid::n_alive() const
{
    int n = 0;
    for(const element& e : elems)
    if(e.alive) ++n;
    return n;
}

double fem_solid::kinetic_energy() const
{
    double ek = 0.0;
    for(int i=0; i<nnode(); ++i)
    ek += 0.5*m[i]*v[i].squaredNorm();
    return ek;
}

double fem_solid::strain_energy() const
{
    double es = 0.0;

    for(int e=0; e<nelem(); ++e)
    {
        const element& el = elems[e];
        if(!el.alive)
        continue;

        const material& mt = mats[el.mat];

        for(int g=0; g<ngp; ++g)
        {
            Mat3 F = Mat3::Zero();
            for(int a=0; a<8; ++a)
            for(int i=0; i<3; ++i)
            for(int J=0; J<3; ++J)
            F(i,J) += x[el.n[a]](i)*(ngp==1 ? dN0[a][J] : dNg[g][a][J]);

            Mat3 E = 0.5*(F.transpose()*F - Mat3::Identity());
            const gpstate& st = gps[e*ngp+g];

            if(mt.type==MAT_J2)
            {
                E(0,0)-=st.Ep[0]; E(1,1)-=st.Ep[1]; E(2,2)-=st.Ep[2];
                E(0,1)-=st.Ep[3]; E(1,0)-=st.Ep[3];
                E(1,2)-=st.Ep[4]; E(2,1)-=st.Ep[4];
                E(2,0)-=st.Ep[5]; E(0,2)-=st.Ep[5];
            }

            const double tr = E.trace();
            Mat3 S = mt.lambda*tr*Mat3::Identity() + 2.0*mt.mu*E;
            es += (1.0-st.d)*0.5*(S.array()*E.array()).sum()*Vel/double(ngp);
        }
    }
    return es;
}

fem_solid::Vec3 fem_solid::support_force() const
{
    return support_force(-1.0e300,1.0e300,-1.0e300,1.0e300,-1.0e300,1.0e300);
}

fem_solid::Vec3 fem_solid::support_force(double x0,double x1,double y0,double y1,double z0,double z1) const
{
    // force of the structure on its supports = load arriving at the fixed dofs
    Vec3 F = Vec3::Zero();
    for(int i=0; i<nnode(); ++i)
    if(fixed[i])
    {
        const Vec3& p = X[i];
        if(p(0)<x0 || p(0)>x1 || p(1)<y0 || p(1)>y1 || p(2)<z0 || p(2)>z1)
        continue;

        Vec3 r = fext[i] + fcon[i] + fcpl[i] + (m[i]-mfl[i])*grav - fint[i];
        for(int d=0; d<3; ++d)
        if(fixed[i] & (1<<d)) F(d) += r(d);
    }
    return F;
}

fem_solid::Vec3 fem_solid::total_load() const
{
    Vec3 F = Vec3::Zero();
    // fluid load: pressure/drag loads, direct-forcing reaction and buoyancy
    for(int i=0; i<nnode(); ++i)
    F += fext[i] + fcpl[i] - mfl[i]*grav;
    return F;
}

double fem_solid::max_vonmises() const
{
    double s = 0.0;
    for(const element& e : elems)
    if(e.alive) s = std::max(s,e.svm);
    return s;
}

double fem_solid::max_displacement() const
{
    double d = 0.0;
    for(int i=0; i<nnode(); ++i)
    d = std::max(d,(x[i]-X[i]).norm());
    return d;
}
